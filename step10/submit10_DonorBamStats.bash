#!/bin/bash

MIN_COVERAGE_DEPTH=1

LEN=$(wc -w < /gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross/donor_list)
FPATH=/gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross_out/1_aligned

for i in $( seq 1 ${LEN} ); do
echo $i
TAG=`awk "NR==$i" /gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross/donor_list`
READ_EXOME=$(samtools view -c -F 4 -L /gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross/targets_chr.bed ${FPATH}/${TAG}_SM_bwa.bam)
READ_SAM=$(samtools view -c -F 4 ${FPATH}/${TAG}_SM_bwa.bam)
COV_EXOME=$(samtools depth -b /gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross/targets_chr.bed ${FPATH}/${TAG}_SM_bwa.bam | wc -l)
COV=$(samtools depth ${FPATH}/${TAG}_SM_bwa.bam | wc -l)
echo $TAG, $READ_EXOME, $READ_SAM, $COV_EXOME, $COV >> ./donor_summary.csv
done
