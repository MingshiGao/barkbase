#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-37%10
#SBATCH --mem=24G
#SBATCH --partition=4hours
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID

mainDir=/data/zusers.ds/gaomingshi/barkbase
chrsize=/data/zusers.ds/gaomingshi/index/canfam3/canfam3.chrom.sizes
sra=$(sed -n ${num}p $mainDir/SRR.list | awk '{print $1}')

sleep $((RANDOM % 300))
mkdir -p /tmp/gaomingshi/${sra}
cd /tmp/gaomingshi/${sra}

samtools view -h $mainDir/nodup/${sra}.nodup.bam | \
  grep -v -e "chrJ" -e "chrA" | \
  samtools view -bS > ${sra}.nodup.bam
sambamba sort -t 10 -n ${sra}.nodup.bam
samtools view ${sra}.nodup.sorted.bam | \
  cut -f 1-11 | sed 'N;s/\n/\t/' | gzip -nc > joined


# Create pseudo replicates bam files
nlines=$(zcat joined | wc -l)
nlines=$(( (nlines + 1) / 2 ))
zcat -f joined | \
  shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(samtools view SRR5427887.nodup.sorted.bam | cut -f 1-11 | wc -c) \
  -nosalt </dev/zero 2>/dev/null) | \
  split -d -l ${nlines} - temp

samtools view -H ${sra}.nodup.sorted.bam > pr1.sam
samtools view -H ${sra}.nodup.sorted.bam > pr2.sam
awk 'BEGIN{OFS="\t"}{
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22
}' temp00 >> pr1.sam
awk 'BEGIN{OFS="\t"}{
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22
}' temp01 >> pr2.sam

# Sort the SAM files and convert them back to BAM format
samtools sort --threads 10 -O BAM pr1.sam > ${sra}.pr1.bam
samtools sort --threads 10 -O BAM pr2.sam > ${sra}.pr2.bam

# Peak calling using MACS2 for both replicates
echo "Calling peaks for $sra..."
macs2 callpeak -t ${sra}.pr1.bam \
  -f BAM -n ${sra}.pr1 \
  -g 2.5e9 -q 0.01 -B --nomodel --nolambda --keep-dup all \
  --call-summits --max-gap 100 --shift -75 --extsize 150

macs2 callpeak -t ${sra}.pr2.bam \
  -f BAM -n ${sra}.pr2 \
  -g 2.5e9 -q 0.01 -B --nomodel --nolambda --keep-dup all \
  --call-summits --max-gap 100 --shift -75 --extsize 150

# Process peaks for both replicates
LC_COLLATE=C sort -k 8gr,8gr ${sra}.pr1_peaks.narrowPeak | \
  awk 'BEGIN{OFS="\t"}{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; print $0}' | \
  gzip -nc > ${sra}.pr1.narrowPeak.gz
LC_COLLATE=C sort -k 8gr,8gr ${sra}.pr2_peaks.narrowPeak | \
  awk 'BEGIN{OFS="\t"}{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; print $0}' | \
  gzip -nc > ${sra}.pr2.narrowPeak.gz


mv ${sra}.pr*.narrowPeak.gz $mainDir/pr/peak/
mv ${sra}.pr*.bam $mainDir/pr/bam

cd $mainDir
rm -rf /tmp/gaomingshi/${sra}
