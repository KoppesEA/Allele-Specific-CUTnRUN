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
DIRsplitbam=./snpsplit/bam
DIRmacs=./MACS/snpsplt_narrowpeaks_wdups_IgG1

#make output directory
mkdir -p $DIRsplitbam

#Load MACS2 2.2.7.1
module load macs/2.2.7.1

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))

#Define snpsplit bam files
fqname=${names[${SLURM_ARRAY_TASK_ID}]}
IP_SPLTBAM_B6=$DIRsplitbam/${fqname}_bt2aln.genome1.bam
IP_SPLTBAM_CAST=$DIRsplitbam/${fqname}_bt2aln.genome2.bam

#Define input IP control could also use [15]
fqname_IgG=${names[7]}
IGG_SPLTBAM_B6=$DIRsplitbam/${fqname_IgG}_bt2aln.genome1.bam
IGG_SPLTBAM_CAST=$DIRsplitbam/${fqname_IgG}_bt2aln.genome1.bam

echo "Starting snpsplit for CnR sample: " $fqname
echo "Control samples is IgG: " $fqname_IgG
echo "IP B6 BAM is: " $IP_SPLTBAM_B6
echo "IgG B6 BAM is: " $IGG_SPLTBAM_B6
echo "IP CAST BAM is: " $IP_SPLTBAM_CAST
echo "IgG B6 BAM is: " $IGG_SPLTBAM_CAST

# Macs2 call with IgG control
macs2 callpeak \
-t $IP_SPLTBAM_B6 \
-c $IGG_SPLTBAM_B6 \
-f BAM \
-g mm \
-q 0.01 \
--bdg \
--keep-dup all \
--outdir $DIRmacs

macs2 callpeak \
-t $IP_SPLTBAM_CAST \
-c $IGG_SPLTBAM_CAST \
-f BAM \
-g mm \
-q 0.01 \
--bdg \
--keep-dup all \
--outdir $DIRmacs

