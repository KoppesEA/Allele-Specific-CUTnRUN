#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRtrim
#SBATCH --output=CnRtrim-%A_%a.txt
#SBATCH --cpus-per-task=4 # 4-cores
#SBATCH --mem=16g # (see also --mem-per-cpu)
#SBATCH --array=0-3 #13 samples

#set path Fastq directories and output directories
DIRfq=./fastq_repcombo
DIRtrim=./trimfastq_repcombo

#make output directory
mkdir $DIRtrim

#Load trimgalore v0.6.5
module load trimgalore/0.6.5

#Extract fq sample names from list of 4 repeated IPs
names=($(printf "IgG\nCTCF\nH3K4me3\nH3K27me3"))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=$DIRfq/BCS_repmerge_${fqname}_1.fastq.gz
INPUT_FASTQ2=$DIRfq/BCS_repmerge_${fqname}_2.fastq.gz

echo "performing quality and adapter trimming for paired-end read on merged replicate R1 and R2:"
echo $INPUT_FASTQ1
echo $INPUT_FASTQ2

#Trim_galore: cutadapt
trim_galore \
--paired \
--illumina \
--quality 20 \
--length 25 \
--e 0.1 \
--phred33 \
--stringency 1 \
--cores 4 \
--gzip \
--fastqc \
--output $DIRtrim \
$INPUT_FASTQ1 $INPUT_FASTQ2

