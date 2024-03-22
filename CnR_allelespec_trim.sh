#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRtrim
#SBATCH --output=CnRtrim-%A_%a.txt
#SBATCH --cpus-per-task=1 # fastqc does not support multithread
#SBATCH --mem=16g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples

#set path Fastq directories and output directories
DIRfq=./rawfastq
DIRtrim=./trimfastq

#make output directory
mkdir $DIRtrim

#Load trimgalore v0.6.5
module load trimgalore/0.6.5

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=${fqname}_1.fastq.gz
INPUT_FASTQ2=${fqname}_2.fastq.gz

echo $INPUT_FASTQ1
echo $INPUT_FASTQ2
echo $OUTPUT

#Trim_galore: cutadapt
mkdir $OUTPUT
trim_galore \
--paired \
--fastqc \
--gzip \
--output $DIRtrim \
$INPUT_FASTQ1 $INPUT_FASTQ2