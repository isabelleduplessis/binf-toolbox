## Get reads mapping length summary stats
```
echo -e "sample\tmin\tq1\tmedian\tq3\tmax\tmean"
for bam in *.bam; do
    echo -e "$(dirname $bam)\t$(python bam_readlen_stats.py $bam | tail -1)"
done
```

## Check strandedness
```
ml RSeQC
infer_experiment.py -i $bamfile -r /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed -s 200000 # sample 200k reads to check
```

## Checking intronic vs exonic for all genes: 

`ml picard`


For hg38
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
```

Get rna seq metrics
```
java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics I=$bamfile O=output_rna_metrics.txt REF_FLAT=refFlat.txt STRAND_SPECIFICITY=NONE # adjust this if you want to use strand specificity
```

Filter for intronic
```
grep -v "^#" output_rna_metrics.txt | \
  awk 'NR==1{for(i=1;i<=NF;i++) col[$i]=i} NR==2{print "PCT_CODING_BASES:", $col["PCT_CODING_BASES"], "\nPCT_INTRONIC_BASES:", $col["PCT_INTRONIC_BASES"]}'
```

## Checking intronic vs exonic for protein coding only:

Get filtered gtfs
```
grep 'gene_type "protein_coding"' /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf | awk '$3=="exon"' > protein_coding_exons.gtf
grep 'gene_type "protein_coding"' /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf | \
  awk '$3=="transcript"' > protein_coding_transcripts.gtf
```


Convert to bed
```
awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' protein_coding_transcripts.gtf > protein_coding_transcripts.bed
awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' protein_coding_exons.gtf > protein_coding_exons.bed
```

Create introns bed
```
bedtools subtract \
  -a protein_coding_transcripts.bed \
  -b protein_coding_exons.bed > protein_coding_introns.bed
```

Get proportions intronic vs exonic (strand specific)
`./count_intron_exon_stranded.sh > intron_exon_results_stranded.tsv`
