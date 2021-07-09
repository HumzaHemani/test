#!/bin/bash

scBAM=${scBAM:default_scBAM} # scBAM for a single cell. 
mutations_csv=${mutations_csv:default_mutations_csv} # reformatted mutations output of step 7 for a single cell. 

while [ $# -gt 0 ]; do            
    if [[ $1 == *"--"* ]]; then
      param="${1/--/}"
      declare $param="$2"
      echo $1 $2 // Optional to see the parameter:value result
    fi       
  shift
done


echo cell: ${scBAM}
echo mutations: ${mutations_csv}

# INITIALIZE VARIABLES/ FUNCTIONS / DIRECTORY PATHS:

extract_meta () {
 	## use awk to find fields that match patterns
 	awk '{ for (i=1; i<=NF; ++i) { if ($i ~ /[DJ]00/) { for (j=1; j<=NF; ++j) { if ($j ~ /CB/) { for (k=1; k<=NF; ++k) { if ($k ~ /UB/) {print $i"___"$j"___"$k} } } } } } }'
}; export -f extract_meta

# GENERATE METADATA:

Rscript /UMI_CORRECTION_4.12.0.R ${mutations_csv}

# EXTRACT METADATA FOR CELL MUTATIONS:

samtools index ${scBAM}

cat ${PWD}/UnfilteredMutations \
| parallel --jobs=30 --max-args=4 samtools view -b -S -h ${scBAM} {1}:{2}-{2} \
'|' java -jar /jvarkit/dist/sam2tsv.jar \
'|' grep -e {4} \
'>>' reads.tsv

cat ${PWD}/UnfilteredMutations \
| parallel --jobs=30 --max-args=4 samtools view ${scBAM} {1}:{2}-{2} \
| extract_meta \
>> meta.tsv

ls -l
