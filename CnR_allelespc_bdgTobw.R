## Erik Koppes
## June 2024
## snpsplt-CutNRun output
## Rscript to convert MACS2 bedgraph output to bigwig
## Some help from Claude (Anthropic AI)

# Load Dependencies 
library(rtracklayer) #[also loads GR, IR and GenomeInfoDb]
library(ChIPpeakAnno) #for toGRanges function [maybe not needed for bdg, but works for peaks]
library(stringr) # for str_remove

## bdg stored with MACS peak calls
MACSbroadpeakIgG1 <- "./MACS/snpsplt_broadpeaks_wdups_IgG1/"
MACSbroadpeakIgG2 <- "./MACS/snpsplt_broadpeaks_wdups_IgG2/"
MACSnarrowpeakIgG1 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG1/"
MACSnarrowpeakIgG2 <- "./MACS/snpsplt_narrowpeaks_wdups_IgG2/"


bdgTobw <- function(bdg_directory) {
  # Import bedGraph files
  for (bedGraph in list.files(bdg_directory, pattern = "treat_pileup\\.bdg$")) {
    print(paste0("Importing: ", bdg_directory, bedGraph))
    bdgGR <- paste0(str_remove(bedGraph, "\\.bdg$"), "_bedGraphGR")
    assign(bdgGR, rtracklayer::import.bedGraph(paste0(bdg_directory, "/", bedGraph),
                                               format = "bedGraph",
                                               genome = Seqinfo(genome="GRCm39")))
  }
  
  # Export to bigwig
  for (bedGraph_name in ls(pattern = "_bedGraphGR$")) {
    print(paste0("Exporting: ", bedGraph_name))
    bedGraph_obj <- get(bedGraph_name)  # Get the actual object
    bwOut <- paste0(bdg_directory, "/", str_remove(bedGraph_name, "_bedGraphGR"), ".bw")
    print(bwOut)
    rtracklayer::export.bw(bedGraph_obj, con = bwOut)
  }
  
  # Clean up: remove bedGraph objects after processing
  rm(list = ls(pattern = "_bedGraph$"), envir = .GlobalEnv)
}

# List of directories
bdg_directories <- c(MACSbroadpeakIgG1, MACSbroadpeakIgG2, MACSnarrowpeakIgG1, MACSnarrowpeakIgG2)

# Process bdg files in each bdg_directories
for (dir in bdg_directories) {
  cat("\nProcessing bdg directory:", dir, "\n")
  bdgTobw(dir)
}