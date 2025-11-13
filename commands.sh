
# Input for workflow would be: .fasta.qz of hg002v1v1.fa.gz, maternal and paternal .gtf

# Split by haplotype: 
grep "MATERNAL" data/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed > data/methylation_maternal.bed
grep "PATERNAL" data/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed > data/methylation_paternal.bed

# Intersect with gene coordinates:
bedtools intersect -a data/maternal.bed -b data/methylation_maternal.bed -wa -wb > data/genes_maternal_overlap.txt
bedtools intersect -a data/paternal.bed -b data/methylation_paternal.bed -wa -wb > data/genes_paternal_overlap.txt

