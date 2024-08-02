# Allele-Specific-CutnRun Analysis on C57B6xCAST mouse ES cells
### Erik Koppes, BIRCWH postdoctoral scholar, Mellissa Mann Lab, Magee Womens Research Institute
### August 2023 - July 2024

## Summary
Code presented here is a series of bash scripts and R Bioconductor analyses to analyze Cut-N-Run chromatin binding assays in an allele-specific manner to determine whether parentally inherited maternal and paternal alleles are bound differently indicative of genomic imprinting. Data and results not shown as this is an unpublished work in progress within the Mann lab.

### Part I SNPsplit B6xCAST genome split
The following commands were run to download the *Mus musculus* Ensembl GRCm39 reference genome and MGP mouse strain variant files to generate *Castaneus* (Cast_EiJ) strain-specific and N-masked genomes. This largely followed the snpsplit genome gen guide (see references) with modification made to run on the University of Pittsburgh HTC.

1. Obtain GRCm39 genome  
  `wget -O Mus_musculus.GRCm39.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz &`

2. Obtain m39 strain specific SNPs in VCF format  
  `wget -O mgp.v5.merged.snps_all.dbSNP142.vcf.gz https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz &`

3. Run SNPsplit_genome_preparation  
   Run the Bash scripts `CastEiJ_snpsplit_genomegen_GRCm39.bash` including output for strain-specific genome

4. Combine chr.txt variant files for CAST  
  From within the output CAST (SNPs_CAST_EiJ) directory run the following interactive commands:   
  `cat chr*.txt > GRCm39_B6CAST_snpsplitSNPs.txt`

### Combine chr.fa files for CAST and N-Masked Reference .fa
From within the output Cast SNPs_CAST_EiJ (CAST_EiJ_full_sequence) and Nmask (CAST_EiJ_N-masked) directories run the following interactive commands:   
`cat chr*.fa > Mus_musculus.GRCm39.109.dna.CASTEiJ.fa`
`cat chr*.fa > Mus_musculus.GRCm39.109.dna.NMASK.fa`

### Part II Cut-N-Run Fastq quality control and pre-processing
1. After manually entering sample basenames as as a single column list in a a plain text text file `CnR_MannLab_fqlist_02262024.txt` run fastqc script `sbatch CnR_allelespec_fastqc.sh` to examine quality of raw fastq data.
2. Summarize pre-trim reports using multiQC on an interactive node:`module load multiqc/1.19` followed by `multiqc --filename "CnR_allele_pretrim_multiQC" ./`
3. Run trimgalore with fastqc output on .fastq input list  `sbatch CnR_allelespec_trim.sh`
4. Summarize post-trim reports using multiQC on an interactive node: `module load multiqc/1.19` followed by `multiqc --filename "CnR_allele_posttrim_multiQC" ./`

## Part III Genomic Alignment
1. Build N-Masked GRCm39 Bowtie2 genomic index: `sbatch CnR_allelespec_bowtie2_NMASK_build.sh`
2. Align Cut-N-Run fastq reads with Bowtie2 to N-Maked reference: `sbatch CnR_allelespec_bowtie2_NMASK.sh`

## Part IV MACS2 Peak Calling
The experiments processed in two sets and utilized two independent IgG. MACS2 gives options for including or discarding "duplicate" reads based on mapped to the same genomic coordianates. Generally for Cut-N-Run using short reads the overlap of mapped reads is great so discarding is not recommended as would reduce quantitative information (as compared to ChIP-Seq which can reduce noise). In addition because reads mapping to the same location may be informative for C57/B6LxCAST SNPs for allele-specifc assays inclusion of duplicate mapped reads with the same coordinates. This is different than optical duplicates that should be removed during QC. MACS2 also allows for either "Narrow peaks" which are generally for TFs and active histone modifications, and for "Broad peaks" that are generally for non-specific DNA binding proteins and repressive histone modifications. Just in case I ran all eight possibilities to compare.

`sbatch CnR_allelespec_MACS2_broadwdupsIgG1con.sh`  
`sbatch CnR_allelespec_MACS2_narrownodupsIgG1con.sh`  
`sbatch CnR_allelespec_MACS2_narrowwdupsIgG1con.sh`  
`sbatch CnR_allelespec_MACS2_broadnodupsIgG1con.sh`  
`sbatch CnR_allelespec_MACS2_broadnodupsIgG2con.sh`  
`sbatch CnR_allelespec_MACS2_broadwdupsIgG2con.sh`  
`sbatch CnR_allelespec_MACS2_narrownodupsIgG2con.sh`  
`sbatch CnR_allelespec_MACS2_narrowwdupsIgG2con.sh`  

## Part V Downstream Bioconductor peak import, annotation, overlap and bedgraph to bigwig conversion
`CnR_allelespec_Rtracklayer_Importbed.R`  
`CnR_allelespec_ChIPPeakAnno.R`  
`CnR_allelespec_ChIPQC.R`  
`CnR_allelespec_ChIPSeakR.R`  
`CnR_allelespc_bdgTobw.R`

## Part VI alternate analyses




## References and Links
https://felixkrueger.github.io/SNPsplit/genome_prep/legacy/
https://ftp.ebi.ac.uk/pub/databases/mousegenomes/
https://www.sanger.ac.uk/data/mouse-genomes-project/
https://www.mousegenomes.org/snps-indels/

