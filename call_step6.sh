#!/bin/bash
# run on R/3.6

# Rscript step6.R BL5481 \
# /home/cifellojp/step5_outs/mutations_NoIntervals \
# /home/cifellojp/step5_in/10x_ref \
# /data/TCR/10X_Genomics/Software/functotator/funcotator_dataSources.v1.6.20190124s \
# ./ \
# ./step6_out


Rscript step6.R BL5481 \
/home/cifellojp/step5/step5_outs/mutations_NoIntervals \
/home/cifellojp/step5/step5_in/10x_ref \
./test_funcotator \
./ \
./step6_out
