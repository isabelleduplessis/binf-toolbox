# Mapping and trimming parameters for ancient nucleic acids

## Bowtie2 vs BWA
bwa aln seems standard for aDNA. No standard yet established for aRNA, but bowtie2 has been successfully used.

## aDNA

Mapping to human genome:
- -l 1024 sets a seed length of 1024, which effectively disables seeding and forces full length alignment fo reads. Increases sensitivity for reads with many mismatches, more accurate alignment if starts of reads are low quality. Slower and more computationally expensive
```
bwa aln -l 1024 -n .01 -o 2
```

Mapping to pathogen genome:
- -n 0.1 allows for more mismatches than -n 0.01
```
bwa aln -n 0.1 -l 1000
```

Competitive mapping to pathogen genomes:
- Stricter than mapping to a single reference genome
```
bwa aln -n 0.01 -l 1000
```

From SPAAM textbook:

Mapping to pathogen (lenient):
```
bwa aln -n 0.01 -l 16
```

Mapping to pathogen (strict):
```
bwa aln -n 0.1 -l 32
```


## aRNA

We see that trimming with a minimum length of 20 can result in poor mapping results for aRNA, to human and bacterial genomes. Use minimum length of 25 or 30. Though for an RNA virus genome we have seen -l 20 work before. You have to check the mapDamage plot after mapping before proceeding.


Mapping to human genome:
- This is what we are working on figuring out

\
Mapping to pathogen genome:
```
bwa aln -n 0.1 -l 32
```


## Mapping QC

MapDamage

## Resources

SPAAM Ancient Metagenomics textbook: https://www.spaam-community.org/intro-to-ancient-metagenomics-book/