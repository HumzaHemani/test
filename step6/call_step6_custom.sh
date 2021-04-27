#!/bin/bash

Rscript /step6.R $1 \
$2 \
$3 \
$4 \
$5 \
$6

for i in $(ls ${5}/annotated/VARIANTS*bash); do
  bash $i
done
