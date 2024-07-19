#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J MACS2spltsnpCnR
#SBATCH --output=MACS2spltsnpCnR-%A_%a.txt
#SBATCH --cpus-per-task=8 # 8-cores per alignment
#SBATCH --mem=24g # (see also --mem-per-cpu)
#SBATCH --array=0-12 #13 samples


#set path Fastq directories and output directories
DIRsam=./bt2align
DIRmacs=./MACS/snpsplt_broadpeaks_nodups_IgG2

#make output directory
mkdir -p $DIRsam

#Load MACS2 2.2.7.1
module load macs/2.2.7.1

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))

#Define snpsplit bam files
fqname=${names[${SLURM_ARRAY_TASK_ID}]}
IP_NMASK_SAM=$DIRsam/${fqname}_bt2aln.genome1.bam


#Define input IP control could also use [6]
fqname_IgG=${names[12]}
IGG_NMASK_SAM_B6=$DIRsam/${fqname_IgG}_bt2aln.genome1.bam
IGG_SPLTBAM_CAST=$DIRsam/${fqname_IgG}_bt2aln.genome2.bam

echo "Starting snpsplit for CnR sample: " $fqname
echo "Control samples is IgG: " $fqname_IgG
echo "IP B6 BAM is: " $IP_NMASK_SAM_B6
echo "IgG B6 BAM is: " $IGG_NMASK_BAM_B6
echo "IP CAST BAM is: " $IP_NMASK_SAM_CAST
echo "IgG B6 BAM is: " $IGG_SPLTBAM_CAST

# Macs2 call with IgG control
macs2 callpeak \
-t $IP_SPLTBAM_B6 \
-c $IGG_NMASK_BAM_B6 \
-f BAMPE \
-g mm \
-q 0.01 \
--bdg \
--broad \
-n ${fqname}_B6_nodupsbroadIgG2con \
--outdir $DIRmacs \


macs2 callpeak \
-t $IP_SPLTBAM_CAST \
-c $IGG_SPLTBAM_CAST \
-f BAMPE \
-g mm \
-q 0.01 \
--bdg \
--broad \
-n ${fqname}_CAST_nodupsbroadIgG2con \
--outdir $DIRmacs

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


