#!/bin/bash

mainDir=/data/zusers/gaomingshi/barkbase
scriptDir=/data/zusers/gaomingshi/barkbase/script

cd $mainDir/Processed-IDRs/
echo "Combining OCRs..."

################################### rOCRs from 24 samples (FRiP > 0.05 & Bladder)
rm -f canFam3-OCRs-manuel.bed
cut -f 1 $mainDir/SRR.filtered.list | while read sra
do
  echo $sra
  cat output.$sra >> canFam3-OCRs-manuel.bed
done

################################## Rescue some Occipital_Cortex and Skeletal Muscle
grep -e Cor -e Skele $mainDir/SRR.list | cut -f 1 | while read sra
do
  awk '$9 > 10' output.$sra >> canFam3-OCRs-manuel.bed
done

input=canFam3-OCRs-manuel.bed
################################## CPKM cutoff - 10th percentile of signal for OCRs
cut -f 5 $input > temp
rm -f temp1
for i in {1..10}
do
  shuf -n 1000000 temp | sort -n | sed -n 100000p >> temp1
done
cutoff=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' temp1)
echo "CPKM cutoff is $cutoff"
rm temp temp1

echo "Filtering OCRs..."
mkdir scratch
awk -v cutoff=$cutoff '$1 !~ /_/ && $3 - $2 >= 150 && $5 >= cutoff && $9 > 5' \
    $input > scratch/tmp.bed
cd scratch

echo "Sorting OCRs..."
sort -k1,1 -k2,2 -V tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

echo "Merging OCRs..."
while [ $num -gt 0 ]
do
  echo -e "\t$num"
  bedtools merge -i sorted -c 4,5 -o collapse,collapse > merge
  python2 $scriptDir/pick.best.peak.py merge > peak-list
  awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
  bedtools intersect -v -a sorted -b rPeaks > remaining
  mv remaining sorted
  num=$(wc -l sorted | awk '{print $1}')
done

mv rPeaks ../tmp.bed
cd ../
rm -r scratch

echo "Accessioning rOCRs..."
sort -k1,1 -k2,2n tmp.bed | awk '{print $0 "\t" "rOCR-"NR}' > rOCRs.bed

