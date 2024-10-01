#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-37
#SBATCH --mem=12G
#SBATCH --cpus-per-task=30
#SBATCH --partition=4hours
#SBATCH --time=4:00:00
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID

mainDir=/data/zusers.ds/gaomingshi/barkbase
ref=/data/zusers.ds/gaomingshi/index/canfam3/bowtie2/canFam3.index

sra=$(sed -n ${num}p $mainDir/SRR.list | awk '{print $1}')
sleep $((RANDOM % 100))

# Trim adapters using Cutadapt
cd $mainDir
mkdir -p fastq_trim
cutadapt --action=mask -m 16 -e 0.2 -j 30 -o fastq_trim/${sra}_1.fastq.gz -p fastq_trim/${sra}_2.fastq.gz fastq/${sra}_1.fastq.gz fastq/${sra}_2.fastq.gz

mkdir -p /tmp/gaomingshi/filter/$sra
cd /tmp/gaomingshi/filter/$sra
bowtie2 --maxins=2000 --threads 30 -x $ref \
        -1 $mainDir/fastq_trim/${sra}_1.fastq.gz \
        -2 $mainDir/fastq_trim/${sra}_2.fastq.gz | \
        samtools view -bS -q 30 -o bam/$sra.bt2.bam
# Filter out unwanted reads
samtools view -F 1804 -f 2 -u $sra.bt2.bam | sambamba sort -o tmp.filt.bam /dev/stdin
# Remove duplicates using Picard
java -Xmx4G -jar /data/zusers.ds/gaomingshi/tools/picard/build/libs/picard.jar MarkDuplicates \
    I=tmp.filt.bam O=tmp.nodup.bam M=$sra.dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
# Remove chrM reads
samtools view -h tmp.nodup.bam | sed '/chrM/d' | samtools view -b -o final.bam /dev/stdin

# Move and index the final BAM file
mv final.bam $sra.nodup.bam
sambamba index $sra.nodup.bam
samtools flagstat $sra.nodup.bam > $sra.nodup.flagstat.qc

# Sort BAM file by name
sambamba sort -t 20 -n $sra.nodup.bam -o ${sra}.nodup.sorted.bam
# Convert BAM to BEDPE format
bedtools bamtobed -bedpe -i ${sra}.nodup.sorted.bam > temp
sort -V -k1,1 -k2,2 temp > ${sra}.bedpe
# Prepare tagAlign format
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' ${sra}.bedpe | \
    awk -F "\t" 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 += 4} else if ($6 == "-") {$3 -= 5} print $0}' > ${sra}.tn5.tagAlign

# Move output files to designated directories
mv ${sra}.tn5.tagAlign $mainDir/tag
mv $sra*hist $mainDir/QC/fragnment_length
mv $sra.nodup.* $mainDir/nodup

cd $mainDir
rm -rf /tmp/gaomingshi/filter/$sra
