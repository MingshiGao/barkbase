#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-37
#SBATCH --mem=12G
#SBATCH --partition=4hours
#SBATCH --time=4:00:00
#SBATCH --output=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.output
#SBATCH --error=/data/zusers.ds/gaomingshi/Logs/jobid_%A_%a.error

num=$SLURM_ARRAY_TASK_ID
mainDir=/data/zusers.ds/gaomingshi/barkbase
sra=$(sed -n ${num}p $mainDir/SRR.list | awk '{print $1}')

singularity exec /home/gaomingshi/encode-chip-seq-pipeline2.simg \
  bash $mainDir/script/run_call_idr.sh $sra

