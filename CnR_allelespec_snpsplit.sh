#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J snpsplt
#SBATCH --output=snpsplt-%A_%a.txt
#SBATCH --cpus-per-task=8 # 8-cores per alignment
#SBATCH --mem=24g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples


#set path Fastq directories and output directories
DIRsam=./bt2align
DIRsplitbam=./snpsplit/bam

#make output directory
mkdir -p $DIRsplitbam

#Load snpsplit v0.6.0
module load snpsplit/0.6.0

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))

#Define Input Trimmed Fastq, Bowtie2 index
fqname=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_SAM=$DIRsam/${fqname}_bt2aln.sam

#Define SNP file
SNPfile=../REF_Sequences/Mus_musculus/Cast_EiJ/CastEiJ_GRCm39/SNPs_CAST_EiJ/GRCm39_B6CAST_snpsplitSNPs.txt

echo "Starting snpsplit for CnR sample: " $fqname
echo "Input sam-file is: " $INPUT_SAM
echo "Input SNP-file is: " $SNPfile

SNPsplit \
--paired \
--output_dir $DIRsplitbam \
--snp_file $SNPfile \
$INPUT_SAM