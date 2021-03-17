#!/bin/bash
# run on R/3.6

Rscript step6.R BL5481 \
~/pipeline/step5/step5_outs/mutations_NoIntervals \
~/pipeline/step5/step5_in/10x_ref \
./funcotator \
./ \
./step6_out

bash ./annotated/VARIANTS_BL5481_1.bash

