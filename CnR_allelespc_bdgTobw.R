## Rscript to convert MACS2 bedgraph output to bigwig

## Additional Script to convert bedgraph to bigwig
## still not working; easier to just make new from .bam

# Load Dependencies 
library(rtracklayer) #[also loads GR, IR and GenomeInfoDb]
library(stringr) # for str_remove

## bdg stored with MACS peak calls
MACSbroadpeakIgG1 <- "./MACS/snpsplt_broadpeaks_wdups_IgG1/"
MACSbroadpeakIgG2 <- "./MACS/snpsplt_broadpeaks_wdups_IgG2/"
MACSnarrowpeakIgG1 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG1/"
MACSnarrowpeakIgG2 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG2/"

# Automated broadPeaks import for BCS allele-specific peaks: Broad-peak-IgG1
for (bedFile in list.files(MACSbroadpeakIgG1, pattern = "treat_pileup.bdg")) {
  assign(paste0(str_remove(bedFile, "treat_pileup.bdg"), ".GR"),
         toGRanges(paste0(MACSbroadpeakIgG1, "/", bedFile), format = "broadPeak"))
  print(paste0(MACSbroadpeakIgG1 ,"/", bedFile))
  print(bedFile)
}
# Automated bed-graph import for Hainer_peaks (from nf-core-MACS pipeline) ##importing cut.bed (filtered output)
MACSdirHainer <- "./Hainer_MACS_oldpipe"
for (bedFile in list.files(MACSdirHainer, pattern = ".clipped.bedGraph")) {
  assign(paste0(str_remove(bedFile, ".clipped.bedGraph"), "_bedGraph"),
         rtracklayer::import.bedGraph(paste0(MACSdirHainer, "/", bedFile)))
  print(paste0(MACSdirHainer ,"/", bedFile))
  print(bedFile)
}

## export to bigwig
for (bedGraph in ls(pattern = "_bedGraph$")) {
  print(bedGraph)
  bwOut<- paste0(MACSdirHainer, "/", str_remove(bedGraph, "_bedGraph"), ".bw")
  print(bwOut)
  rtracklayer::export.bw(bedGraph, con = bwOut, format = "bw", genome="mm10")
}

## import+export function
for (bedGraph in list.files(MACSdirHainer, pattern = ".clipped.bedGraph")) {
  print(bedGraph)
  bwOut<- paste0(MACSdirHainer, "/", str_remove(bedGraph, ".clipped.bedGraph"), ".bw")
  print(bwOut)
  rtracklayer::wigToBigWig(x = paste0(MACSdirHainer, "/", bedGraph),
                           dest = bwOut, seqinfo= SeqInfoMm)
}

library(BSgenome)
available.genomes()
BiocManager::install("BSgenome.BSgenome.Mmusculus.UCSC.mm10") #install as needed
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Sscrofa.UCSC.susScr11)
Sscr11_BSgenome <- BSgenome.Sscrofa.UCSC.susScr11
seqnames(Sscr11_BSgenome)

# Specify the path to genome FASTA file
Mm10_fa_file <- "./Ben_MACS_nfcorv321/igv/genome.fa"

# Import Sscr11 .fa (Ensembl) reference pig genome v11) using Biostrings::readDNAStringSet()
Mm10_fa <- readDNAStringSet(Mm10_fa_file)
SeqInfoMm <- seqinfo(Mm10_fa)
seqinfo(IgG_R1_bedGraph)

Sscr11_BSgenome <- BSgenome.Sscrofa.UCSC.susScr11
seqnames(Sscr11_BSgenome)
test <- seqlengths((Sscr11_BSgenome))
