#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J BT2_GRCm39
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem=128g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=BT2_GRCm39_indexing.out

#set path nmask fa file and bt2 index output directory
REFDIR=../REF_Sequences/Mus_musculus/GRCm39_ref
GENOME=$REFDIR/Mus_musculus.GRCm39.dna.toplevel.fa
OUTPUT=$REFDIR/BT2index
 
module load gcc/8.2.0
module load bowtie2/2.3.4.2 

bowtie2-build \
--threads 16 \
$GENOME \
$OUTPUT


