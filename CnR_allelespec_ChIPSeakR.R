## Annotation with ChIP Seekr
## See Online: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html
## Install ChIPseeker, tidygraph and TxDb.Mm if needed
# BiocManager::install("ChIPseeker") 
# BiocManager::install("tidygraph") #had to update for ChIPSeeker
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")

## Load Dependencies
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#library(EnsDb.Hsapiens.v75)
#library(clusterProfiler)
#library(AnnotationDbi)
#library(org.Hs.eg.db) Mus.musculus

# Set txdb
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
keytypes(TxDb.Mmusculus.UCSC.mm10.ensGene)
#[1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"     "TXNAME"  

## Sample MetaData
BCS_meta<- data.frame(
  Sample_ID = paste0("BS", seq(1:6)),
  Sample_Name =   c("NUP107_ab85916", "NUP107_ab236634", "CTCF", "H3k27me3", "H3k4me3", "IgG"))

# Load Ben nf-core Pipeline v.3.2.1 data
v321igvdirBen <- "./Ben_MACS_nfcorv321/igv"
samplefilesBCS <- list.files(v321igvdirBen, pattern = ".relaxed.bed", full.names=T)
names(samplefilesBCS) <- BCS_meta$Sample_Name[1:5] #no IgG peak file output


# Annotate BCS peaks
peakAnnoListBCS <- lapply(samplefilesBCS, annotatePeak, TxDb=txdb, 
                                              tssRegion=c(-1000, 1000), verbose=FALSE)
# Graph BCS peaks
peakAnnoListBCS
plotAnnoBar(peakAnnoListBCS)
plotDistToTSS(peakAnnoListBCS, title="Distribution of ChIP Peaks relative to TSS")

## Load Hainer data and annotate features
MACSdirHainer <- "./Hainer_MACS_oldpipe"
samplefilesHainer <- list.files(MACSdirHainer, pattern = ".narrowPeak", full.names=T)
names(samplefilesHainer) <- basename(samplefilesHainer)

peakAnnoListHainer <- lapply(samplefilesHainer, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoListHainer
plotAnnoBar(peakAnnoListHainer)
plotDistToTSS(peakAnnoListHainer, title="Distribution of ChIP Peaks relative to TSS")





# ## Additional Plots
# # See https://hbctraining.github.io/Intro-to-ChIPseq/lessons/chipseeker_visualization.html
# # Assign peak data to variables
# BSC_H3K27me3 <- samplefilesBCS[4]
# BSC_H3K4me3 <- samplefilesBCS[5]
# 
# # Plot covplot
# covplot(BSC_H3K27me3, weightCol="V5")
# covplot(BSC_H3K4me3, weightCol="V5")
# 
# ## Plot promoter TSS
# # Prepare the promotor regions
# promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
# # Calculate the tag matrix
# tagMatrixList <- lapply(as.list(c(BSC_H3K27me3, BSC_H3K4me3)), getTagMatrix, windows=promoter) ##would use full BSC_ but empyt for Nup107
# ## Profile plots
# plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")
# 
# # Plot heatmap
# tagHeatmap(tagMatrixList)
#            #, xlim=c(-1000, 1000), color=NULL)
