##  Rscript ChIP-Peak analysis for Nup107 peaks
##  ChIPQC, 
##  Erik Koppes
##  May 05 2024

library(ChIPQC)
## remove ChIPpeakAnno for ggplot2 unload--
detach("package:ChIPpeakAnno", unload=TRUE)
library(ChIPQC)


## From::https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/ChIP-seq/doc/ChIP-Seq_step.by.step.html#/33

peaks<-ls(pattern=".GR")

ChIP_Samples <- 
  data.frame(
    Sample_ID = c(paste0("BCS", LETTERS[1:6]), paste0("BCS", seq(1:7))),
    CellType = "mESC",
    Factor = c("NUP107_4", "NUP107_6", "NUP107_8", "NUP107_10", "H3K4me", "IgG",
               "CTCF_rep1", "CTCF_rep2", "H3K4me3_rep1", "H3K4me3_rep2", "H3K27me3_rep1", "H3K27me3_rep2", "IgG"),
    bamReads = 
    Peaks = peaks
  )
,
Sample_Name =   c("NUP107_4", "NUP107_6", "NUP107_8", "NUP107_10", "H3K4me", "IgG",
                  "CTCF_rep1", "CTCF_rep2", "H3K4me3_rep1", "H3K4me3_rep2", "H3K27me3_rep1", "H3K27me3_rep2", "IgG"))
samples <- 
  data.frame(SampleID=gsub("^(.*?).q1.*$", "\\1", names(peaks)), 
             Factor=gsub("^(.*?)_.*$", "\\1", names(peaks)),
             Replicate=gsub("^.*?_rep(.*?).q1.*$", "\\1", 
                            names(peaks)),
             bamReads=gsub("^(.*?).q1.*$", "\\1.rmdup.bam", 
                           names(peaks)),
             Peaks=names(peaks),
             PeakCaller="narrow")
samples[1:3, ]

## before filter
exampleExp <- ChIPQC(experiment=samples, annotaiton="hg19")
ChIPQCreport(exampleExp, reportFolder="ChIPQC", facetBy="Factor")
## after IDR filter
samples.fil <- samples
samples.fil$Peaks <- gsub("^(.*?).q1.*$", "\\1.fil.bed", names(peaks))
samples.fil$PeakCaller="bed"
exampleExp.fil <- ChIPQC(experiment=samples.fil, annotaiton="hg19")
ChIPQCreport(exampleExp.fil, 
             reportFolder="ChIPQC.fil", 
             facetBy="Factor")
## filtered by qValue from the MACS2 result
samples.fil.q01 <- samples.fil
samples.fil.q01$Peaks <- 
  gsub("^(.*?).q1.*$", "\\1.fil.q01.bed", names(peaks))
exampleExp.fil.q01 <- ChIPQC(experiment=samples.fil.q01, 
                             annotaiton="hg19")
