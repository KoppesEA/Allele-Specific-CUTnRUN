#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 6-00:00 # 6d Max Runtime in D-HH:MM
#SBATCH -J CnR_BBinterleave_dddup_unpairepair
#SBATCH --output=CnR_BBMap_clumpify_dedup_unpairepair-%A_%a.txt
#SBATCH --cpus-per-task=1 # 1-core, not parallel
#SBATCH --mem=8g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples

#set path Fastq directories and output directories
DIRfq=./rawfastq
DIRdedup=./fastq_clumpify_dedup_unpairrepair

#make output directory
mkdir $DIRdedup

#Load BBMAP v38.51
module load bbmap/38.51

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1=$DIRfq/${fqname}_1.fastq.gz
INPUT_FASTQ2=$DIRfq/${fqname}_2.fastq.gz

OUTPUT_FASTQ1=$DIRdedup/${fqname}_optdedupe_unpairrepair_1.fastq.gz
OUTPUT_FASTQ2=$DIRdedup/${fqname}_optdedupe_unpairrepair_2.fastq.gz

echo "performing optical read deduplication with BBMap clumpify dedupe"
echo $INPUT_FASTQ1
echo $INPUT_FASTQ2
echo $OUTPUT_FASTQ1
echo $OUTPUT_FASTQ2


# CnR sequencing on NextSeq at Pitt Core

clumpify.sh \
in1=$INPUT_FASTQ1 \
in2=$INPUT_FASTQ2 \
out1=$OUTPUT_FASTQ1 \
out2=$OUTPUT_FASTQ2 \
groups=auto \
dedupe=t \
optical=t \
dupedist=40 \
spany=t \
adjacent=t \
unpair=t \
repair=t

