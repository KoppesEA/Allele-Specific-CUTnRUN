##  Rscript ChIP-Peak analysis for Nup107 peaks
##  Erik Koppes
##  May 2 2024

#install BioCManager if Needed; install ChIPpeakAnno and biomaRt; open librariers
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("biomaRt")

library(ChIPpeakAnno)
library(biomaRt)


## feature annotation
listMarts()
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                dataset="mmusculus_gene_ensembl")
Annotation <- getAnnotation(mart,
                            featureType = c("TSS", "miRNA", "Exon", "5utr", "3utr", "ExonPlusUtr", "transcript"))
peak_list <- c(BS1_R1.seacr.peaks.relaxed_bedImport, BS2_R1.seacr.peaks.relaxed_bedImport, BS3ctcf_R1.seacr.peaks.relaxed_bedImport, BS4k27_R1.seacr.peaks.relaxed_bedImport, BS5k4_R1.seacr.peaks.relaxed_bedImport,
               )
peak1_annot <- 
  annotatePeakInBatch(
  CTCFrep2_R1.seacr.peaks.stringent_bedImport,
  mart,
  featureType = c("TSS"),
  AnnotationData = Annotation,
  output = c("nearestLocation", "overlapping", "both", "shortestDistance", "inside",
             "upstream&inside", "inside&downstream", "upstream", "downstream",
             "upstreamORdownstream", "nearestBiDirectionalPromoters"),
  multiple = c(TRUE, FALSE),
  maxgap = -1L,
  PeakLocForDistance = c("start", "middle", "end", "endMinusStart"),
  FeatureLocForDistance = c("TSS", "middle", "start", "end", "geneEnd"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  bindingRegion = NULL
)

y = peak1_annot$distancetoFeature[!is.na(peak1_annot$distancetoFeature)&annotatedPeak$fromOverlappingOrNearest == "NearestStart"]
hist(y, xlab = "Distance To Nearest TSS", main = "", breaks = 1000, xlim = c(min(y)-100, max(y) + 100))

peak1_annot <- 
  annotatePeakInBatch(
    H3k27merep2_R1.seacr.peaks.stringent_bedImport,
    mart,
    featureType = c("TSS"),
    AnnotationData = Annotation,
    output = c("nearestLocation", "overlapping", "both", "shortestDistance", "inside",
               "upstream&inside", "inside&downstream", "upstream", "downstream",
               "upstreamORdownstream", "nearestBiDirectionalPromoters"),
    multiple = c(TRUE, FALSE),
    maxgap = -1L,
    PeakLocForDistance = c("start", "middle", "end", "endMinusStart"),
    FeatureLocForDistance = c("TSS", "middle", "start", "end", "geneEnd"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = TRUE,
    bindingRegion = NULL
  )

y = peak1_annot$distancetoFeature[!is.na(peak1_annot$distancetoFeature)]
hist(y, xlab = "Distance To Nearest TSS", main = "", breaks = 1000, xlim = c(min(y)-100, max(y) + 100))

peak1_annot <- 
  annotatePeakInBatch(
    H3k4me_R1.seacr.peaks.stringent_bedImport,
    mart,
    featureType = c("TSS"),
    AnnotationData = Annotation,
    output = c("nearestLocation", "overlapping", "both", "shortestDistance", "inside",
               "upstream&inside", "inside&downstream", "upstream", "downstream",
               "upstreamORdownstream", "nearestBiDirectionalPromoters"),
    multiple = c(TRUE, FALSE),
    maxgap = -1L,
    PeakLocForDistance = c("start", "middle", "end", "endMinusStart"),
    FeatureLocForDistance = c("TSS", "middle", "start", "end", "geneEnd"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = TRUE,
    bindingRegion = NULL
  )

y = peak1_annot$distancetoFeature[!is.na(peak1_annot$distancetoFeature)]
hist(y, xlab = "Distance To Nearest TSS", main = "", breaks = 1000, xlim = c(min(y)-100, max(y) + 100))
