#!/bin/bash

sra=$1

chrsize=/data/zusers.ds/gaomingshi/index/canfam3/canfam3.chrom.sizes
mainDir=/data/zusers/gaomingshi/barkbase
pr1=$mainDir/pr/peak/${sra}.pr1.narrowPeak.gz
pr2=$mainDir/pr/peak/${sra}.pr2.narrowPeak.gz
peak=$mainDir/peak/narrowpeak/${sra}.narrowPeak.gz
tag=$mainDir/tag/${sra}.tn5.tagAlign

mkdir -p /tmp/gaomingshi/$sra
cd /tmp/gaomingshi/$sra

python3 $(which encode_task_idr.py) \
        $pr1 $pr2 $peak \
        --prefix $sra \
        --idr-thresh 0.05 \
        --peak-type narrowPeak \
        --idr-rank signal.value \
        --fraglen 73 \
        --chrsz $chrsize \
        --regex-bfilt-peak-chr-name 'chr[\dXY]+' \
        --ta $tag

mv $sra.idr0.05.narrowPeak.gz $mainDir/idr
mv $sra.idr0.05.bfilt.frip.qc $mainDir/idr

cd $mainDir
rm -r /tmp/gaomingshi/$sra
