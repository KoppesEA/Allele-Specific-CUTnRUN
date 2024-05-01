#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRfqc
#SBATCH --output=CnR_BBdedup_fqc-%A_%a.txt
#SBATCH --cpus-per-task=1 # fastqc does not support multithread
#SBATCH --mem=16g # prob overkill (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples

#set path Fastq directories and output directories
DIRfq=./fastq_clumpify_dedup
DIRfqQC=$DIRfq/FastQC

#make output directory
mkdir $DIRfqQC

#Load fastqc v0.11.9
module load fastqc/0.11.9

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=${fqname}_optdedupe_1.fastq.gz
INPUT_FASTQ2=${fqname}_optdedupe_2.fastq.gz
echo $INPUT_FASTQ1
echo $INPUT_FASTQ2

#perform a FastQC analysis for every fastq, including for each read pair
fastqc -o $DIRfqQC $DIRfq/$INPUT_FASTQ1
fastqc -o $DIRfqQC $DIRfq/$INPUT_FASTQ2
