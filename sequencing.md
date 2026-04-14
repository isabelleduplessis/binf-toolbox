PCR duplicates: Same DNA is amplified resulting in multiple identical reads
Optical duplciates: One read is sequenced multiple times


original:  M-8-1.bt2.sort.bam
57266780


java -jar $EBROOTPICARD/picard.jar MarkDuplicates     I=M-8-1.bt2.sort.bam     O=seqs0.bam REMOVE_SEQUENCING_DUPLICATES=true     M=marked_dup_metrics.txt
57266780


java -jar $EBROOTPICARD/picard.jar MarkDuplicates     I=M-8-1.bt2.sort.bam     O=seqs_dups0.bam REMOVE_DUPLICATES=true     M=marked_dup_metrics.txt
1514849




"The need to de-duplicate datasets, in which identical reads are collapsed to a single read, is much debated amongst RNA researchers because of uncertainty about whether duplicates represent biological expression or an artefact of the PCR process [26]. Considering the short nature of our RNA reads and the generally high duplication rate, we surmised that these were more likely to be a PCR artefact than a reflection of biological expression. We performed all analyses using both the unmodified (‘duplicated’) and de-duplicated sets and found that in all cases, de-duplicated data made more biological sense."




original:  120916_R.bt2.sort.bam
2183364


java -jar $EBROOTPICARD/picard.jar MarkDuplicates     I=120916_R.bt2.sort.bam     O=seqs0.bam REMOVE_SEQUENCING_DUPLICATES=true     M=marked_dup_metrics.txt
2183364


java -jar $EBROOTPICARD/picard.jar MarkDuplicates     I=120916_R.bt2.sort.bam     O=seqs_dups0.bam REMOVE_DUPLICATES=true     M=marked_dup_metrics.txt
227591