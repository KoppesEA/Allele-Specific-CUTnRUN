#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRtrim
#SBATCH --output=CnRtrim-%A_%a.txt
#SBATCH --cpus-per-task=4 # 4-cores
#SBATCH --mem=16g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples

#set path Fastq directories and output directories
DIRfq=./rawfastq
DIRtrim=./trimfastq

#make output directory
mkdir $DIRtrim

#Load Bowtie2 v2.4.5
module load bowtie2/2.4.5

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=${fqname}_1.fastq.gz
INPUT_FASTQ2=${fqname}_2.fastq.gz

echo "performing quality and adapter trimming for paired-end reads:"
echo $INPUT_FASTQ1
echo $INPUT_FASTQ2

# Bowtie2 with CUT&RUN parameters
bowtie2 \
-p 4 \
-x $BOWTIE_INDEX \
-U $INPUT_FASTQ_CUTADAPT \
-S ./Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}.sam