#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-37
#SBATCH --mem=12G
#SBATCH --partition=4hours
#SBATCH --time=04:00:00
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID

mainDir=/data/zusers.ds/gaomingshi/barkbase
chrsize=/data/zusers.ds/gaomingshi/index/canfam3/canfam3.chrom.sizes
file=$mainDir/SRR.list

sra=$(sed -n "${num}p" $file | awk '{print $1}')

sleep $((RANDOM % 100))

mkdir -p /tmp/gaomingshi/${sra}
cd /tmp/gaomingshi/${sra}
cp $mainDir/nodup/${sra}.nodup.bam .

# Use only autosomal reads for peak calling with MACS2
samtools view -h ${sra}.nodup.bam | grep -v -e "chrJ" -e "chrA" | samtools view -bS > ${sra}.nodup.filtered.bam
macs2 callpeak -t ${sra}.nodup.filtered.bam \
    -f BAM -n ${sra} \
    -g 2.5e9 -q 0.01 -B --nomodel --nolambda --keep-dup all \
    --call-summits --max-gap 100 \
    --shift -75 --extsize 150 --SPMR

# Process peak files
LC_COLLATE=C sort -k 8gr,8gr ${sra}_peaks.narrowPeak | \
    awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; if ($2<0) $2=0; if ($3<0) $3=0; print $0}' | \
    gzip -nc > ${sra}.narrowPeak.gz
zcat ${sra}.narrowPeak.gz | bedClip stdin $chrsize temp
sort -V -k1,1 -k2,2 temp | bedtools merge -i stdin -c 8,9 -o max,max | \
    awk -v OFS="\t" '{print $1,$2,$3,"'${sra}'_peak_"NR,$4,$5 }' > ${sra}.uniq.peak.bed

# Generate files for track hub
LC_COLLATE=C sort -k1,1 -k2,2n ${sra}.uniq.peak.bed | \
    awk -v OFS="\t" '{if ($6 > 2) print $1,$2,$3,$4}' > temp1
bedToBigBed temp1 $chrsize ${sra}.bb
LC_COLLATE=C sort -k1,1 -k2,2n ${sra}_treat_pileup.bdg > sorted.bdg
bedGraphToBigWig sorted.bdg $chrsize ${sra}_treat_pileup.bw

mv ${sra}.uniq.peak.bed $mainDir/peak
mv ${sra}.narrowPeak.gz *xls $mainDir/peak/narrowpeak
mv *.bb $mainDir/hub
mv ${sra}.*bw $mainDir/hub

cd $mainDir
rm -rf /tmp/gaomingshi/${sra}
