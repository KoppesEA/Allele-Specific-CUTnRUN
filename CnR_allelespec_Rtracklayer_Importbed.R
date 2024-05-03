## Rscript to import .bedgraph files and compare Heiner et al C&R with new BSC data
## C&R nf-core pipe Seacr peakcalling outputs a bedgraph like file
## https://github.com/FredHutch/SEACR
## Erik Koppes 2/8/2024


# Load Dependencies 
library(rtracklayer) #[also loads GR, IR and GenomeInfoDb]
library(stringr)
library(ChIPpeakAnno) # for EZ Venndigrams
library(ChIPQC)
BiocManager::install("ChIPQC")
## Sample MetaData
BCS_meta<- data.frame(
  Sample_ID = c(paste0("BCS", LETTERS[1:6]), paste0("BCS", seq(1:7))),
  Sample_Name =   c("NUP107_4", "NUP107_6", "NUP107_8", "NUP107_10", "H3K4me", "IgG",
                    "CTCF_rep1", "CTCF_rep2", "H3K4me3_rep1", "H3K4me3_rep2", "H3K27me3_rep1", "H3K27me3_rep2", "IgG"))


MACSbroadpeakIgG1 <- "./MACS/snpsplt_broadpeaks_wdups_IgG1/"
MACSbroadpeakIgG2 <- "./MACS/snpsplt_broadpeaks_wdups_IgG2/"

# Automated broadPeaks import for BCS allele-specific peaks: Broad-peak-IgG1
for (bedFile in list.files(MACSbroadpeakIgG1, pattern = ".broadPeak")) {
  assign(paste0(str_remove(bedFile, ".broadPeak"), ".GR"),
  toGRanges(paste0(MACSbroadpeakIgG1, "/", bedFile), format = "broadPeak"))
  print(paste0(MACSbroadpeakIgG1 ,"/", bedFile))
  print(bedFile)
}

# Automated broadPeaks import for BCS allele-specific peaks: Broad-peak-IgG2
for (bedFile in list.files(MACSbroadpeakIgG2, pattern = ".broadPeak")) {
  assign(paste0(str_remove(bedFile, ".broadPeak"), ".GR"),
         toGRanges(paste0(MACSbroadpeakIgG2, "/", bedFile), format = "broadPeak"))
  print(paste0(MACSbroadpeakIgG2 ,"/", bedFile))
  print(bedFile)
}

# Automated broadPeaks import for BCS allele-specific peaks: Broad-peak-IgG1
for (bedFile in list.files(MACSnarrowpeakIgG1, pattern = ".narrowPeak")) {
  assign(paste0(str_remove(bedFile, ".narrowPeak"), ".GR"),
         toGRanges(paste0(MACSnarrowpeakIgG1, "/", bedFile), format = "narrowPeak"))
  print(paste0(MACSnarrowpeakIgG1 ,"/", bedFile))
  print(bedFile)
}

# Automated broadPeaks import for BCS allele-specific peaks: Broad-peak-IgG2
for (bedFile in list.files(MACSnarrowpeakIgG2, pattern = ".narrowPeak")) {
  assign(paste0(str_remove(bedFile, ".narrowPeak"), ".GR"),
         toGRanges(paste0(MACSnarrowpeakIgG2, "/", bedFile), format = "narrowPeak"))
  print(paste0(MACSnarrowpeakIgG2 ,"/", bedFile))
  print(bedFile)
}

## Examine mat/pat CTCF overlap for rep1 and rep2
peak_overlap_CTCF_rep1 <- findOverlaps(BCS1_B6_wdupsbroadIgG1con_peaks.GR, BCS1_CAST_wdupsbroadIgG1con_peaks.GR)

## Use makeVennDiagram Function from ChipPeakAnno to generate Venn Diagram 

## Should use narrowPeak and IgG2 for CTCF and H3K4me3
# Venndiagram of CTCF BCS rep1 and 2 B6(mat) and CAST(pat)
makeVennDiagram(list(BCS1_B6_wdupsbroadIgG1con_peaks.GR, BCS1_CAST_wdupsbroadIgG1con_peaks.GR,BCS2_B6_wdupsbroadIgG1con_peaks.GR, BCS2_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("CTCF_rep1_B6", "CTCF_rep1_CAST", "CTCF_rep2_B6", "CTCF_rep2_CAST"),
                totalTest=181125,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))


# Venndiagram of H3K4me3 BCS rep1 and 2 B6(mat) and CAST(pat)
makeVennDiagram(list(BCS3_B6_wdupsbroadIgG1con_peaks.GR, BCS3_CAST_wdupsbroadIgG1con_peaks.GR,BCS4_B6_wdupsbroadIgG1con_peaks.GR, BCS4_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("H3K4me3_rep1_B6", "H3K4me3_rep1_CAST", "H3K4me3_rep2_B6", "H3K4me3_rep2_CAST"),
                totalTest=102832,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of H3K27me3 BCS rep1 and 2 B6(mat) and CAST(pat)
makeVennDiagram(list(BCS5_B6_wdupsbroadIgG1con_peaks.GR, BCS5_CAST_wdupsbroadIgG1con_peaks.GR,BCS6_B6_wdupsbroadIgG1con_peaks.GR, BCS6_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("H3K27me3_rep1_B6", "H3K27me3_rep1_CAST", "H3K27me3_rep2_B6", "H3K27me3_rep2_CAST"),
                totalTest=102832,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))


# Comparison of Nup107 antibodies
# maximum levels for makeVenndiagram =5; so do 4 f-way comparisons
# Venndiagram of Nup107 BCS _4 and _6 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSA_B6_wdupsbroadIgG1con_peaks.GR, BCSA_CAST_wdupsbroadIgG1con_peaks.GR,BCSB_B6_wdupsbroadIgG1con_peaks.GR, BCSB_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_4_B6", "NUP107_4_CAST", "NUP107_6_B6", "NUP107_6_CAST"),
                totalTest=127919,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of Nup107 BCS _4 and _8 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSA_B6_wdupsbroadIgG1con_peaks.GR, BCSA_CAST_wdupsbroadIgG1con_peaks.GR,BCSC_B6_wdupsbroadIgG1con_peaks.GR, BCSC_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_4_B6", "NUP107_4_CAST", "NUP107_8_B6", "NUP107_8_CAST"),
                totalTest=97933,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of Nup107 BCS _4 and _10 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSA_B6_wdupsbroadIgG1con_peaks.GR, BCSA_CAST_wdupsbroadIgG1con_peaks.GR,BCSD_B6_wdupsbroadIgG1con_peaks.GR, BCSD_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_4_B6", "NUP107_4_CAST", "NUP107_10_B6", "NUP107_10_CAST"),
                totalTest=102832,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of Nup107 BCS _6 and _8 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSB_B6_wdupsbroadIgG1con_peaks.GR, BCSB_CAST_wdupsbroadIgG1con_peaks.GR,BCSC_B6_wdupsbroadIgG1con_peaks.GR, BCSC_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_6_B6", "NUP107_6_CAST", "NUP107_8_B6", "NUP107_8_CAST"),
                totalTest=67644,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of Nup107 BCS _6 and _10 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSB_B6_wdupsbroadIgG1con_peaks.GR, BCSB_CAST_wdupsbroadIgG1con_peaks.GR,BCSD_B6_wdupsbroadIgG1con_peaks.GR, BCSD_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_6_B6", "NUP107_6_CAST", "NUP107_10_B6", "NUP107_10_CAST"),
                totalTest=64914,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))

# Venndiagram of Nup107 BCS _8 and _10 antibody B6(mat) and CAST(pat)
makeVennDiagram(list(BCSC_B6_wdupsbroadIgG1con_peaks.GR, BCSC_CAST_wdupsbroadIgG1con_peaks.GR,BCSD_B6_wdupsbroadIgG1con_peaks.GR, BCSD_CAST_wdupsbroadIgG1con_peaks.GR),
                NameOfPeaks=c("NUP107_8_B6", "NUP107_8_CAST", "NUP107_10_B6", "NUP107_10_CAST"),
                totalTest=64914,scaled=FALSE, euler.d=FALSE,  ## set totalTest= total of 2 reps; tested with 250k then calc
                fill=c("#FFBAD2", "#66CCFF", "#F48FB1","#0099CC"), # circle fill color
                col=c("#F20056", "#003399","#F20056", "#003399"), #circle border color
                cat.col=c("#FFBAD2", "#66CCFF", "#F48FB1", "#0099CC"))
