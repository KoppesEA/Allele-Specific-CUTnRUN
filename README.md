# Allele-Specific-CutnRun Analysis on C57B6xCAST mouse ES cells
### Erik Koppes, BIRCWH postdoctoral scholar, Mellissa Mann Lab, Magee Womens Research Institute

## Part I snp-split B6xCAST genome split

## Summarry
The following commands were run to download the GRCm39 genomes and then to generate Castenous (Cast_EiJ) strain-specific genomes. While this largely followed the snpsplit genome gen guide (see references) there were some specific issues getting the backdated GRCm38 to work. This was run July 26th 2023 on the PITT HTC.

### Obtain GRCm39 genome
`wget -O Mus_musculus.GRCm39.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz &`

### Obtain m39 strain specific SNPs in VCF format
`wget -O mgp.v5.merged.snps_all.dbSNP142.vcf.gz https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz &`

### Run SNPsplit_genome_preparation
Run the Bash scripts `CastEiJ_snpsplit_genomegen_GRCm39.bash` including output for strain-specific genome

### Combine chr.txt variant files for CAST
From within the output CAST (SNPs_CAST_EiJ) directory run the following interactive commands:   
`cat chr*.txt > GRCm39_B6CAST_snpsplitSNPs.txt`

### Combine chr.fa files for CAST and N-Masked Reference .fa
From within the output Cast SNPs_CAST_EiJ (CAST_EiJ_full_sequence) and Nmask (CAST_EiJ_N-masked) directories run the following interactive commands:   
`cat chr*.fa > Mus_musculus.GRCm39.109.dna.CASTEiJ.fa`
`cat chr*.fa > Mus_musculus.GRCm39.109.dna.NMASK.fa`


## References and Links
https://felixkrueger.github.io/SNPsplit/genome_prep/legacy/
https://ftp.ebi.ac.uk/pub/databases/mousegenomes/
https://www.sanger.ac.uk/data/mouse-genomes-project/
https://www.mousegenomes.org/snps-indels/

