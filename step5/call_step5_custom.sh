Rscript step5.R \
$1 \
$2 \
$3 \
$4 \
$5 \
$6 \
$7

for i in $(ls ${6}/VARIANTS*bash); do
  bash $i
done
