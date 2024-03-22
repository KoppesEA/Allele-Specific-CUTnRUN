#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRtrim
#SBATCH --output=CnRtrim-%A_%a.txt
#SBATCH --cpus-per-task=8 # 8-cores per alignment
#SBATCH --mem=32g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples

#set path Fastq directories and output directories
DIRtrim=./trimfastq
DIRbam=./bt2align

#make output directory
mkdir $DIRtrim

#Load Bowtie2 v2.4.5
module load bowtie2/2.4.5

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

#Define Input Trimmed Fastq, Bowtie2 index
fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=$DIRtrim/${fqname}_1.fastq.gz
INPUT_FASTQ2=$DIRtrim/${fqname}_2.fastq.gz
BOWTIE_INDEX=../REF_Sequences/Mus_musculus/Cast_EiJ/CastEiJ_GRCm39/
SAMout=$DIRbam/${fqname}_bt2aln.sam

echo "performing bowtie2 readmapping to GRCm39 N-masked C57/B6N x CastEiJ SNPs"
echo $INPUT_FASTQ1
echo $INPUT_FASTQ2
echo $READ_BASE
echo $INPUT_FASTQ_CUTADAPT
echo $BOWTIE_INDEX

# Bowtie2 with CUT&RUN parameters
bowtie2 \
-p 8 \
-x $BOWTIE_INDEX \
-U $INPUT_FASTQ_CUTADAPT \
--end-to-end \
--very-sensitive \
--no-mixed \
--no-discordant \
-I 10 -X 700 \
--dovetail \
-1 INPUT_FASTQ1 \
-2 INPUT_FASTQ2 \
-S $SAMout



-S ./Bowtie2cutUMIAdaptalignedm18_v99/${READ_BASE}.sam