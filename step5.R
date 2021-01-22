
version <- '4.6.3'

# Load Libraries
library(tidyverse)
library(R.utils)
library(argparse)
# source('./parameters')

## use these when testing
# args <- c()
# args$decom_meta <- "./decom_meta_1000cells.csv"
# args$aligned_dir <- "./aligned"
# args$scBAM_dir <- "./scBAM"
# args$snp_out <- "./snps"
# args$reference <- "/data/TCR/__REFERENCE_GENOMES/10X_2020A"
# args$scripts_dir <- "."
# args$out <- "/data/TCR/10X_Genomics/data/test/singularity_jeffrey"

parser <- ArgumentParser()
parser$add_argument("decom_meta", help="csv containing columns: barcode,BC_ID,Person,Age,Sex,BL_ID,Code,visit,study,hash")
parser$add_argument("aligned_dir", help="Directory containing alignment files.")
parser$add_argument("scBAM_dir", help="Directory containing BAM files for each cell.")
parser$add_argument("snp_out", help="SNP results directory from previous step.")
parser$add_argument("reference", help="10x reference directory.")
parser$add_argument("scripts_dir", help="Where to write the bash files calling gatk.")
parser$add_argument("out", help="Where all pipeline outs will be saved.")
args <- parser$parse_args()

# process passed args
decom_meta_hash <- read.csv(args$decom_meta)
scripts3 <- file.path(args$scripts_dir, 'variant_calling', '4.6.3_GenerateScripts_VariantCalling_ReadBased_NoIntervals')
scripts_dir <- args$scripts_dir
aligned <- args$aligned_dir
scBAM <- args$scBAM_dir
out <- args$out
ref <- args$reference

var_out <- file.path(out,'mutations_NoIntervals')
ase_out <- file.path(out, 'allele_specific')
snp_out <- args$snp_out

# Load Data
decom_meta_hash$X1 <- paste(decom_meta_hash$barcode, '1', sep = '-')

master_list <- c()
dir.create(file.path(scripts_dir, 'variant_calling'))
dir.create(scripts3)
tag_dir <- file.path(scripts3,'TL')
dir.create(tag_dir)
dir.create(var_out)
samples <- unique(decom_meta_hash$BL_ID)
for (i in samples) dir.create(file.path(var_out,i))

# Generate Scripts for Variant-Calling
for (sam in samples) {
  sample_cells <- decom_meta_hash[decom_meta_hash$BL_ID == sam,]$X1
  spl <- split(sample_cells, ceiling(seq_along(sample_cells)/215))
  for (ls in c(1:length(spl))) {
    write(unlist(spl[ls]), file = paste(tag_dir, '/TL_', sam, '_', ls, sep = ''))
    
    bash <- paste('#!/bin/bash

# LOAD MODULES:

module load parallel
module load bamtools
module load samtools
module load GATK/4.1.0.0

# SAMPLE NUMBER:
# TAG LIST NUMBER:

SAMPLE=',sam,'
TL=',ls,'
VERSION=',version,'

# INITIALIZE VARIABLES/ DIRECTORY PATHS:
cd ', getwd(), '

TNAME=CB
DIR=',out,'/2_scBAM
scBAM=', scBAM,'/',sam,'
DIR_SCRIPT=',script_dir,'/variant_calling/4.6.3_GenerateScripts_VariantCalling_ReadBased_NoIntervals
DIR_O=',var_out,'
BIGBAM=',aligned,'/${SAMPLE}_SM_bwa.bam
BIGBAM_LOCAL=/lscratch/${SLURM_JOB_ID}/${SAMPLE}_SM_bwa.bam
REF=',ref,'
REF_LOCAL=/lscratch/${SLURM_JOB_ID}/genome.fa
SNPs=',snp_out,'
REF_VAR_LOCAL=/lscratch/${SLURM_JOB_ID}/',sam,'_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa
TMPDIR=/lscratch/${SLURM_JOB_ID}/

# COPY BAM/REFERENCE FILE TO LOCAL SCRATCH:

cp $BIGBAM /lscratch/${SLURM_JOB_ID}/. || exit 1
cp ${BIGBAM}.bai /lscratch/${SLURM_JOB_ID}/. || exit 1
cp $REF/genome.fa /lscratch/${SLURM_JOB_ID}/. || exit 1
cp $REF/genome.dict /lscratch/${SLURM_JOB_ID}/. || exit 1
cp $REF/genome.fa.fai /lscratch/${SLURM_JOB_ID}/. || exit 1
cp $SNPs/',sam,'_SM_bwa_RawSNPs_FLTR_SNP_consensus.fa /lscratch/${SLURM_JOB_ID}/. || exit 1

# READ IN TAGS - PIPE INTO GNU PARALLEL :
#						SPLIT BAM INTO SINGLE CELLS > CALL VARIANTS > FILTER/AGGREGATE VARIANTS to mutations.csv

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --progress --jobs 20 cp ${scBAM}/${SAMPLE}_{}.bam /lscratch/${SLURM_JOB_ID}/.

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --progress --jobs 20 samtools index \\
/lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}.bam
 
cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --progress --jobs 20 gatk --java-options "\'-Xmx1G\'" AddOrReplaceReadGroups \\
-I /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}.bam \\
-O /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_UMI_SM.bam \\
-ID ${SAMPLE}_{} \\
-LB MissingLibrary \\
-PL ILLUMINA \\
-PU ${SAMPLE}_{} \\
-SM ${SAMPLE}_{}

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --progress --jobs 20 samtools index \\
/lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_UMI_SM.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs 12 gatk --java-options "\'-Xmx8g -XX:+UseConcMarkSweepGC\'" SplitNCigarReads \\
-R ${REF_LOCAL} \\
-I /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_UMI_SM.bam \\
-O /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_UMI_SM_ST.bam

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs 12 gatk --java-options "\'-Xmx8g -XX:+UseConcMarkSweepGC\'" Mutect2 \\
-R ${REF_VAR_LOCAL} \\
-I /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_UMI_SM_ST.bam \\
-I ${BIGBAM_LOCAL} \\
-tumor ${SAMPLE}_{} \\
-normal ${SAMPLE}_combined \\
-DF MappingQualityAvailableReadFilter \\
-DF MappingQualityReadFilter \\
-DF MappingQualityNotZeroReadFilter \\
-O /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_var.vcf

cat ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL} \\
| parallel --jobs 16 gatk --java-options "\'-Xmx4g -XX:+UseConcMarkSweepGC\'" FilterMutectCalls \\
-V /lscratch/${SLURM_JOB_ID}/${SAMPLE}_{}_var.vcf \\
-O ${DIR_O}/${SAMPLE}/${SAMPLE}_{}_var_FLTR.vcf \\
--tumor-lod 5.3 \\
--disable-tool-default-read-filters \\
--min-median-base-quality 0 \\
--min-median-mapping-quality 0 \\
--min-median-read-position 0

LEN=$(wc -w < ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL})

for i in $( seq 1 ${LEN} ); do
TAG=`awk "NR==$i" ${DIR_SCRIPT}/TL/TL_${SAMPLE}_${TL}`
VART=$(cat /lscratch/${SLURM_JOB_ID}/${SAMPLE}_${TAG}_var.vcf | wc -l)
VAR=$(cat  ${DIR_O}/${SAMPLE}/${SAMPLE}_${TAG}_var_FLTR.vcf | grep "PASS" | wc -l)
READ_SAM=$(samtools view -c -F 4 /lscratch/${SLURM_JOB_ID}/${SAMPLE}_${TAG}.bam)
UMI=$(samtools view  /lscratch/${SLURM_JOB_ID}/${SAMPLE}_${TAG}.bam | grep -o \'UB:............\' | grep -o \'..........$\' | uniq -c | wc -l)
COV=$(samtools mpileup /lscratch/${SLURM_JOB_ID}/${SAMPLE}_${TAG}.bam | awk -v X="${MIN_COVERAGE_DEPTH}" \'$4>=X\' | wc -l)
echo $SAMPLE, $TAG, $VART, $VAR, $READ_SAM, $UMI, $COV >> ${DIR_O}/mutations.csv
done', sep = '')
    
    write(bash, file = paste(scripts3, '/VARIANTS_', sam, '_', ls, '.bash', sep = ''))
    
    master_list <- c(master_list, paste(sam,'_',ls, sep = ''))
  }
}

ml_swarm <- paste('bash ./VARIANTS_',master_list,'.bash', sep='')
write(ml_swarm, file = paste(scripts3,'/swarm.swarm',sep=''))

swarm_bash <- paste(
'#!/bin/bash
#swarm -f swarm.swarm --job-name=mut_exome --time=36:00:00 --partition=norm -g 4 -t 2 --gres=lscratch:10 --logdir=./__logs
swarm -f swarm.swarm --job-name=mut_exome --time=68:00:00 --partition=norm -g 30 -t 18 --gres=lscratch:110 --logdir=./__logs'
,sep='')

write(swarm_bash, file = paste(scripts3,'/swarm.bash',sep=''))
