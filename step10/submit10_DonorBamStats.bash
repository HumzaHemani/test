#!/bin/bash

donor_list=${sample:def_donor_list}
exome_sc=${exome_sc:def_exome_sc}
output_summary=${output_summary:def_output_summary}

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 #// Optional to see the parameter:value result
    fi       
  shift
done

LEN=$(wc -w < ${donor_list})

for i in $( seq 1 ${LEN} ); do
echo $i
TAG=`awk "NR==$i" /gpfs/gsfs5/users/TCR/10X_Genomics/code/test/mutations_cross/donor_list`
READ_EXOME=$(samtools view -c -F 4 -L /targets_chr.bed ${exome_sc}/${TAG}_SM_bwa.bam)
READ_SAM=$(samtools view -c -F 4 ${exome_sc}/${TAG}_SM_bwa.bam)
COV_EXOME=$(samtools depth -b /targets_chr.bed ${exome_sc}/${TAG}_SM_bwa.bam | wc -l)
COV=$(samtools depth ${exome_sc}/${TAG}_SM_bwa.bam | wc -l)
echo $TAG, $READ_EXOME, $READ_SAM, $COV_EXOME, $COV >> ${output_summary}
done
