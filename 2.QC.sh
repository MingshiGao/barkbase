#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-37
#SBATCH --mem=10G
#SBATCH --partition=4hours
#SBATCH --time=04:00:00
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID

mainDir=/data/zusers.ds/gaomingshi/barkbase
tss=/data/zusers.ds/gaomingshi/index/canfam3/4000.10bp.fantom.cage.tss.bed
sra=$(sed -n ${num}p $mainDir/SRR.list | awk '{print $1}')

mkdir -p /tmp/gaomingshi/$sra
cd /tmp/gaomingshi/$sra
sleep $((RANDOM % 100))

### examine TSS-enrichment
cp $mainDir/hub/${sra}_treat_pileup.bw .
bigWigAverageOverBed ${sra}_treat_pileup.bw $tss ${sra}.out
python3 $mainDir/script/tss_coverage_bw.py ${sra}.out 400 > $mainDir/QC/tss-enrich/${sra}.tss.summary

### Calculate the FRiP score
cp $mainDir/tag/${sra}.tn5.tagAlign .
cut -f 1-3 $mainDir/peak/${sra}.uniq.peak.bed | uniq | bedtools merge -i stdin > ${sra}.temp.bed
bedtools intersect -wa -u -a ${sra}.tn5.tagAlign -b ${sra}.temp.bed > temp
n1=$(wc -l ${sra}.tn5.tagAlign | awk '{print $1}')
n2=$(wc -l temp | awk '{print $1}')
frip=$(echo "$n1 $n2" | awk '{printf "%.5f", $2 / $1}')
echo -e "$sra\t$n1\t$n2\t$frip" >> $mainDir/QC/frip_autosome.txt

cd $mainDir
rm -rf /tmp/gaomingshi/$sra

