#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J BT2_BxCnmask
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem=128g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=BT2_BxCnmask_indexing.out

#set path nmask fa file and bt2 index output directory
BxCDIR=../REF_Sequences/Mus_musculus/Cast_EiJ/CastEiJ_GRCm39/CAST_EiJ_N-masked
GENOME=$BxCDIR/Mus_musculus.GRCm39.109.dna.NMASK.fa
OUTPUT=$BxCDIR/BT2index
 
module load gcc/8.2.0
module load bowtie2/2.3.4.2 

bowtie2-build \
--threads 16 \
$GENOME \
$OUTPUT


