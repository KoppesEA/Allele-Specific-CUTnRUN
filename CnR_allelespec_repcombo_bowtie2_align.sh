#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnRaln
#SBATCH --output=CnRalnCombo-%A_%a.txt
#SBATCH --cpus-per-task=8 # 8-cores per alignment
#SBATCH --mem=32g # (see also --mem-per-cpu)
#SBATCH --array=0-3 #13 samples

#set path Fastq directories and output directories
DIRtrim=./trimfastq_repcombo
DIRbam=./bt2align_repcombo

#make output directory
mkdir $DIRbam

#Load Bowtie2 v2.4.5
module load bowtie2/2.4.5

#Extract fq sample names from list text file
names=($(printf "IgG\nCTCF\nH3K4me3\nH3K27me3"))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

#Define Input Trimmed Fastq, Bowtie2 index
fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=$DIRtrim/BCS_repmerge_${fqname}_1_val_1.fq.gz
INPUT_FASTQ2=$DIRtrim/BCS_repmerge_${fqname}_2_val_2.fq.gz
BOWTIE_INDEX=../REF_Sequences/Mus_musculus/Cast_EiJ/CastEiJ_GRCm39/CAST_EiJ_N-masked/BT2index
SAMout=$DIRbam/BCS_repmerge_${fqname}_bt2aln.sam

echo "performing bowtie2 readmapping to GRCm39 N-masked C57/B6N x CastEiJ SNPs"
echo "fastq CnR sample: " $fqname
echo "fastq trimmed read1 file: "$INPUT_FASTQ1
echo "fastq trimmed read2 file: "$INPUT_FASTQ2
echo "Bowtie2 Index: " $BOWTIE_INDEX

# Bowtie2 with CUT&RUN parameters
bowtie2 \
-p 8 \
-x $BOWTIE_INDEX \
--end-to-end \
--very-sensitive \
--no-mixed \
--no-discordant \
-I 10 -X 700 \
--dovetail \
-1 $INPUT_FASTQ1 \
-2 $INPUT_FASTQ2 \
-S $SAMout