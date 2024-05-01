#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRfqc
#SBATCH --output=CnRfqc-%A_%a.txt
#SBATCH --cpus-per-task=1 # fastqc does not support multithread
#SBATCH --mem=16g # prob overkill (see also --mem-per-cpu)
#SBATCH --array=0-3 #13 samples

#set path Fastq directories and output directories
DIRfq=./fastq_repcombo
DIRfqQC=$DIRfq/FastQC

#make output directory
mkdir $DIRfqQC

#Load fastqc v0.11.9
module load fastqc/0.11.9

#Extract fq sample names from list of 4 repeated IPs
names=($(printf "IgG\nCTCF\nH3K4me3\nH3K27me3"))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=BCS_repmerge_${fqname}_1.fastq.gz
INPUT_FASTQ2=BCS_repmerge_${fqname}_2.fastq.gz

#perform a FastQC analysis for every fastq, including for each read pair
fastqc -o $DIRfqQC $DIRfq/$INPUT_FASTQ1
fastqc -o $DIRfqQC $DIRfq/$INPUT_FASTQ2