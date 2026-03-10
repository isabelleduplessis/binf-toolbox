#!/bin/bash

# Below uses strand speceficity, but remove -s to not check that

INTRONS=/fh/fast/blancomelo_d/user/iduplessis/paleolung/ID010_aRNA_hmap_optz/featureCounts/protein_coding_introns.bed
EXONS=/fh/fast/blancomelo_d/user/iduplessis/paleolung/ID010_aRNA_hmap_optz/featureCounts/protein_coding_exons.bed

# Print header
echo -e "sample\tintronic_reads\texonic_reads\ttotal\tfraction_intronic\tfraction_exonic"

for bam in /fh/fast/blancomelo_d/user/iduplessis/paleolung/ID010_aRNA_hmap_optz/samples/*001/bwa.sorted.dedup.q30.bam; do
        sample=$(basename $(dirname $bam))
    
    INTRONIC=$(bedtools intersect \
        -a $bam \
        -b $INTRONS \
        -u -s | samtools view -c)
    
    EXONIC=$(bedtools intersect \
        -a $bam \
        -b $EXONS \
        -u -s | \
    bedtools intersect \
        -a stdin \
        -b $INTRONS \
        -v -s | samtools view -c)
    
    TOTAL=$((INTRONIC + EXONIC))
    
    FRAC_INTRONIC=$(echo "scale=3; $INTRONIC / $TOTAL" | bc)
    FRAC_EXONIC=$(echo "scale=3; $EXONIC / $TOTAL" | bc)
    
    echo -e "${sample}\t${INTRONIC}\t${EXONIC}\t${TOTAL}\t${FRAC_INTRONIC}\t${FRAC_EXONIC}"

done
