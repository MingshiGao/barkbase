#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-20
#SBATCH --mem=36G
#SBATCH --partition=4hours
#SBATCH --time=04:00:00
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID

mainDir=/data/zusers.ds/gaomingshi/barkbase
pubDir=/data/public_html_users/gaomingshi/barkbase2021/canfam3
chrsize=/data/zusers.ds/gaomingshi/index/canfam3/canFam3.chrom.sizes

sra=$(sed -n ${num}p $mainDir/SRR.list | awk '{print $1}')

sleep $((RANDOM % 150))
mkdir -p /tmp/gaomingshi/$sra
cd /tmp/gaomingshi/$sra

### Generate cutsite BED and RPM signal BigWig ###
cp $mainDir/nodup/$sra.nodup.bam .
bedtools bamtobed -i $sra.nodup.bam | grep "/1" | \
  awk -v OFS="\t" '$1 !~ /chrJ|chrA/ {if ($6 == "+") print $1,$2+4,$2+5; else if ($6 == "-") print $1,$3-6,$3-5}' \
  > $sra.cutsite.bed
bedtools genomecov -bga -scale $n -i $sra.cutsite.bed -g $chrsize > temp.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n temp.bedGraph > sort.bedGraph
bedGraphToBigWig sort.bedGraph $chrsize $sra.RPM.bw

### Calculate CPKM signal and z-score ###
zcat $mainDir/idr/$sra.idr0.05.narrowPeak.gz | cut -f 1-3 | sort -k1,1 -k2,2 -Vu | awk -v sra=$sra '{OFS="\t"}{print $1,$2,$3,sra"_"NR}' > bed
bigWigAverageOverBed $sra.RPM.bw bed out.tab
python2 $mainDir/script/calculate.zscore.sh out.tab | sort -k1,1 > zscore
paste bed zscore | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8 "\t" "." "\t" $7 "\t" 1 "\t" $5}' > output.$sra

mv $sra.RPM.bw $pubDir/$sra.cutsite.RPM.bw
mv $sra.cutsite.bed $mainDir/cutsite
mv output.$sra $mainDir/Processed-IDRs/

cd $mainDir
rm -rf /tmp/gaomingshi/$sra

