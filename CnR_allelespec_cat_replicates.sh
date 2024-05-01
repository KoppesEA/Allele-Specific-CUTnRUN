#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J CnR_repCat
#SBATCH --output=CnR_repCat-%A_%a.txt
#SBATCH --cpus-per-task=8 #
#SBATCH --mem=16g # (see also --mem-per-cpu)

DIRfq=./rawfastq
DIRfq_combo=./fastq_repcombo

#make output directory
mkdir $DIRfq_combo

#Extract fq sample names from list text file
names=($(cat CnR_MannLab_fqlist_02262024.txt))

#Define fq names
fqname_IgG_1=${names[5]}
fqname_IgG_2=${names[12]}
fqname_CTCF_1=${names[6]}
fqname_CTCF_2=${names[7]}
fqname_H3K4me3_1=${names[8]}
fqname_H3K4me3_2=${names[9]}
fqname_H3K27me3_1=${names[10]}
fqname_H3K27me3_2=${names[11]}
echo "CnR Sample replicate names: "
echo $fqname_IgG_1 $fqname_IgG_2
echo $fqname_CTCF_1 $fqname_CTCF_2
echo $fqname_H3K4me3_1 $fqname_H3K4me3_2
echo $fqname_H3K27me3_1 $fqname_H3K27me3_2

#Define input fq r1 and r2 files
IgG_1_FASTQ1=$DIRfq/${fqname_IgG_1}_1.fastq.gz
IgG_1_FASTQ2=$DIRfq/${fqname_IgG_1}_2.fastq.gz
IgG_2_FASTQ1=$DIRfq/${fqname_IgG_2}_1.fastq.gz
IgG_2_FASTQ2=$DIRfq/${fqname_IgG_2}_2.fastq.gz
CTCF_1_FASTQ1=$DIRfq/${fqname_CTCF_1}_1.fastq.gz
CTCF_1_FASTQ2=$DIRfq/${fqname_CTCF_1}_2.fastq.gz
CTCF_2_FASTQ1=$DIRfq/${fqname_CTCF_2}_1.fastq.gz
CTCF_2_FASTQ2=$DIRfq/${fqname_CTCF_2}_2.fastq.gz
H3K4me3_1_FASTQ1=$DIRfq/${fqname_H3K4me3_1}_1.fastq.gz
H3K4me3_1_FASTQ2=$DIRfq/${fqname_H3K4me3_1}_2.fastq.gz
H3K4me3_2_FASTQ1=$DIRfq/${fqname_H3K4me3_2}_1.fastq.gz
H3K4me3_2_FASTQ2=$DIRfq/${fqname_H3K4me3_2}_2.fastq.gz
H3K27me3_1_FASTQ1=$DIRfq/${fqname_H3K27me3_1}_1.fastq.gz
H3K27me3_1_FASTQ2=$DIRfq/${fqname_H3K27me3_1}_2.fastq.gz
H3K27me3_2_FASTQ1=$DIRfq/${fqname_H3K27me3_2}_1.fastq.gz
H3K27me3_2_FASTQ2=$DIRfq/${fqname_H3K27me3_2}_2.fastq.gz
echo "Merging Fastq Replicate Paired-end Reads: "
echo $IgG_1_FASTQ1 $IgG_1_FASTQ2
echo $IgG_2_FASTQ1 $IgG_2_FASTQ2
echo $CTCF_1_FASTQ1 $CTCF_1_FASTQ2
echo $CTCF_2_FASTQ1 $CTCF_2_FASTQ2
echo $H3K4me3_1_FASTQ1 $H3K4me3_1_FASTQ2
echo $H3K4me3_2_FASTQ1 $H3K4me3_2_FASTQ2
echo $H3K27me3_1_FASTQ1 $H3K27me3_1_FASTQ2
echo $H3K27me3_2_FASTQ1 $H3K27me3_2_FASTQ2

#Concatenate gzipped replicate FASTQ1 and FASTQ2 read files
cat $IgG_1_FASTQ1 $IgG_2_FASTQ1 > $DIRfq_combo/BCS_repmerge_IgG_1.fastq.gz
cat $IgG_1_FASTQ2 $IgG_2_FASTQ2 > $DIRfq_combo/BCS_repmerge_IgG_2.fastq.gz
cat $CTCF_1_FASTQ1 $CTCF_2_FASTQ1 > $DIRfq_combo/BCS_repmerge_CTCF_1.fastq.gz
cat $CTCF_1_FASTQ2 $CTCF_2_FASTQ2 > $DIRfq_combo/BCS_repmerge_CTCF_2.fastq.gz
cat $H3K4me3_1_FASTQ1 $H3K4me3_2_FASTQ1 > $DIRfq_combo/BCS_repmerge_H3K4me3_1.fastq.gz
cat $H3K4me3_1_FASTQ2 $H3K4me3_2_FASTQ2 > $DIRfq_combo/BCS_repmerge_H3K4me3_2.fastq.gz
cat $H3K27me3_1_FASTQ1 $H3K27me3_2_FASTQ1 > $DIRfq_combo/BCS_repmerge_H3K27me3_1.fastq.gz
cat $H3K27me3_1_FASTQ2 $H3K27me3_2_FASTQ2 > $DIRfq_combo/BCS_repmerge_H3K27me3_2.fastq.gz