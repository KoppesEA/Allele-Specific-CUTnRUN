## Erik Koppes
## June 2024
## snpsplt-CutNRun output
## Rscript to convert MACS2 bedgraph output to bigwig

# Load Dependencies 
library(rtracklayer) #[also loads GR, IR and GenomeInfoDb]
library(ChIPpeakAnno) #for toGRanges function [maybe not needed for bdg, but works for peaks]
library(stringr) # for str_remove

## bdg stored with MACS peak calls
MACSbroadpeakIgG1 <- "./MACS/snpsplt_broadpeaks_wdups_IgG1/"
MACSbroadpeakIgG2 <- "./MACS/snpsplt_broadpeaks_wdups_IgG2/"
MACSnarrowpeakIgG1 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG1/"
MACSnarrowpeakIgG2 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG2/"

# Automated bed-graph import for snpsplit-MACS peak calls ##treat_pileup.bdg (Treated Sample normalized to IgG), with designation of GRCm39 annotation
for (bedGraph in list.files(MACSbroadpeakIgG1, pattern = "treat_pileup.bdg")) {
  print(paste0(MACSbroadpeakIgG1, bedGraph))
  print(bedGraph)
  bdgGR <- paste0(str_remove(bedGraph, ".clipped.bedGraph"), "_bedGraph")
  # print(bdgGR)
  assign(bdgGR, rtracklayer::import.bedGraph(paste0(MACSbroadpeakIgG1, "/", bedGraph),
                                             format = "bedGraph",
                                             genome = Seqinfo(genome="GRCm39")))
}


## export to bigwig
for (bedGraph_name in ls(pattern = "_bedGraph$")) {
  print(bedGraph_name)
  bedGraph_obj <- get(bedGraph_name)  # Get the actual object
  bwOut <- paste0(MACSbroadpeakIgG1, "/", str_remove(bedGraph_name, "_bedGraph"), ".bw")
  print(bwOut)
  rtracklayer::export.bw(bedGraph_obj, con = bwOut)  # No need to specify format="bw" here
}

# 
# 
# for (bedGraph in ls(pattern = "_bedGraph$")) {
#   print(bedGraph)
#   bwOut<- paste0(MACSbroadpeakIgG1, "/", str_remove(bedGraph, "_bedGraph"), ".bw")
#   print(bwOut)
#   rtracklayer::export.bw(bedGraph, con = bwOut, format = "bw") ## genome="mm10"
# }
# 
# # BCS4_B6_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph@seqnames@values <- paste0("chr", BCS4_B6_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph@seqnames@values) dont't set names, but set seqinfo
# rtracklayer::export.bw(BCS4_B6_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph, con = "./BCS4_B6_wdupsbroadIgG1con_treat_pileup.bw", format = "bw")
# rtracklayer::export.bw(BCS1_B6_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph, con = "./BCS1_B6_wdupsbroadIgG1con_treat_pileup.bw", format = "bw")
# 
# BCS4_B6_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph@seqinfo = Seqinfo(genome="GRCm39")
# 
# BCSF_CAST_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph
# BCSF_CAST_wdupsbroadIgG1con_treat_pileup.bdg_bedGraph
# ####################-----------#####################

# ## import+export function
# for (bedGraph in list.files(MACSdirHainer, pattern = ".clipped.bedGraph")) {
#   print(bedGraph)
#   bwOut<- paste0(MACSdirHainer, "/", str_remove(bedGraph, ".clipped.bedGraph"), ".bw")
#   print(bwOut)
#   rtracklayer::wigToBigWig(x = paste0(MACSdirHainer, "/", bedGraph),
#                            dest = bwOut, seqinfo= SeqInfoMm)
# }
# 
# library(BSgenome)
# available.genomes()
# BiocManager::install("BSgenome.BSgenome.Mmusculus.UCSC.mm10") #install as needed
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(BSgenome.Sscrofa.UCSC.susScr11)
# Sscr11_BSgenome <- BSgenome.Sscrofa.UCSC.susScr11
# seqnames(Sscr11_BSgenome)
# 
# # Specify the path to genome FASTA file
# Mm10_fa_file <- "./Ben_MACS_nfcorv321/igv/genome.fa"
# 
# # Import Sscr11 .fa (Ensembl) reference pig genome v11) using Biostrings::readDNAStringSet()
# Mm10_fa <- readDNAStringSet(Mm10_fa_file)
# SeqInfoMm <- seqinfo(Mm10_fa)
# seqinfo(IgG_R1_bedGraph)
# 
# Sscr11_BSgenome <- BSgenome.Sscrofa.UCSC.susScr11
# seqnames(Sscr11_BSgenome)
# test <- seqlengths((Sscr11_BSgenome))
