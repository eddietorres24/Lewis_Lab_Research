### Script for ChIPQC

#Sample code for troubleshooting
##show contrasts for a DBA
#dba.show(DBA, bContrast=T)

## Load Libraries
library(BiocManager)
library(ChIPQC)
library(ggplot2)
library(DiffBind)
library(tidyverse)

## Load Sample Data
setwd("C:\\Users\\eddie\\Research\\GitHub")
K27_abc_samples <- read.csv("./DiffBind_CAF-1_K27_abc_ET.csv")
K27_CS_samples <- read.csv("./DiffBind_CAF-1_K27_CS_ET.csv")
K36_samples <- read.csv("./DiffBind_CAF-1_K36_ET.csv")
H4K20_samples <- read.csv("./DiffBind_CAF-1_H4K20_ET.csv")
H3K4me2_samples <- read.csv("./DiffBind_CAF-1_H3K4me2_ET.csv")
ATAC_samples <- read.csv("./DiffBind_CAF-1_ATAC_ET.csv")

## Make ChIPQC Object
register(SerialParam())
chipObj <- ChIPQC(K27_CS_samples) 

## Make QC Report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_H3K27me3_CS", reportFolder="ChIPQCreport")


### DiffBind (scratch)

#Read in sample sheet containing CAF-1 mutant data *If running ChIP QC, you can use "samples" variable
caf_dba_K27_abc <- dba(sampleSheet = K27_abc_samples)
caf_dba_K27_CS <- dba(sampleSheet = K27_CS_samples)
caf_dba_K36 <- dba(sampleSheet = K36_samples)
caf_dba_H4K20 <- dba(sampleSheet = H4K20_samples)
caf_dba_H3K4me2 <- dba(sampleSheet = H3K4me2_samples)
caf_dba_ATAC <- dba(sampleSheet = ATAC_samples)

##Plot correlation plot between all samples included on data sheet, save in a file
CorrPlot_K27_abc <- plot(caf_dba_K27_abc)
CorrPlot_K27_CS <- plot(caf_dba_K27_CS)
CorrPlot_K36 <- plot(caf_dba_K36)
CorrPlot_H4K20 <- plot(caf_dba_H4K20)
CorrPlot_H3K4me2 <- plot(caf_dba_H3K4me2)
CorrPlot_ATAC <- plot(caf_dba_ATAC)

#Using RPKM
CorrPlot_RPKM <- dba.plotHeatmap(caf_dba_K27_CS, score=DBA_SCORE_RPKM_FOLD)

#Count Reads
caf_dba_K27_abc <- dba.count(caf_dba_K27_abc)
caf_dba_K27_CS <- dba.count(caf_dba_K27_CS)
caf_dba_K36 <- dba.count(caf_dba_K36)
caf_dba_H4K20 <- dba.count(caf_dba_H4K20)
caf_dba_H3K4me2 <- dba.count(caf_dba_H3K4me2)
caf_dba_ATAC <- dba.count(caf_dba_ATAC, summits=75)

##Plot new correlation plot based on count data, save in a file
CorrPlot_count_K27_abc <- plot(caf_dba_K27_abc)
CorrPlot_count_K27_CS <- plot(caf_dba_K27_CS)
CorrPlot_count_K36 <- plot(caf_dba_K36)
CorrPlot_count_H4K20 <- plot(caf_dba_H4K20)
CorrPlot_count_H3K4me2 <- plot(caf_dba_H3K4me2)
CorrPlot_count_ATAC <- plot(caf_dba_ATAC)

#Normalization
caf_dba_K27_abc_norm <- dba.normalize(caf_dba_K27_abc, normalize="RLE")
caf_dba_K27_CS_norm <- dba.normalize(caf_dba_K27_CS, normalize="lib")
caf_dba_K36_norm <- dba.normalize(caf_dba_K36, normalize="RLE")
caf_dba_H4K20_norm <- dba.normalize(caf_dba_H4K20, normalize="RLE")
caf_dba_H3K4me2_norm <- dba.normalize(caf_dba_H3K4me2, normalize="RLE")
caf_dba_ATAC_norm <- dba.normalize(caf_dba_ATAC, normalize=DBA_NORM_NATIVE)

#Model design & Contrast (what comparisons do you want to make?)
caf_dba_K27_abc_norm <- dba.contrast(caf_dba_K27_abc_norm, categories=DBA_FACTOR, minMembers = 2)
caf_dba_K27_CS_norm <- dba.contrast(caf_dba_K27_CS_norm, categories=DBA_FACTOR, minMembers = 2)
caf_dba_K36_norm <- dba.contrast(caf_dba_K36_norm, categories=DBA_FACTOR, minMembers = 2)
caf_dba_H4K20_norm <- dba.contrast(caf_dba_H4K20_norm, categories=DBA_FACTOR, minMembers = 2)
caf_dba_H3K4me2_norm <- dba.contrast(caf_dba_H3K4me2_norm, categories=DBA_FACTOR, minMembers = 2)
caf_dba_ATAC_norm <- dba.contrast(caf_dba_ATAC_norm, categories=DBA_FACTOR, minMembers = 2)

#Blacklist
#caf_dba <- dba.blacklist(caf_dba, blacklist=FALSE, greylist=FALSE)

#Differential Analysis
caf_dba_K27_abc_norm <- dba.analyze(caf_dba_K27_abc_norm)
caf_dba_K27_CS_norm <- dba.analyze(caf_dba_K27_CS_norm)
caf_dba_K36_norm <- dba.analyze(caf_dba_K36_norm)
caf_dba_H4K20_norm <- dba.analyze(caf_dba_H4K20_norm)
caf_dba_H3K4me2_norm <- dba.analyze(caf_dba_H3K4me2_norm)
caf_dba_ATAC_norm <- dba.analyze(caf_dba_ATAC_norm)

#Plot PCA for All mods (normalized peak calls)
dba.plotPCA(caf_dba_K27_abc_norm, label=DBA_FACTOR)
dba.plotPCA(caf_dba_K27_CS_norm, label=DBA_FACTOR)
dba.plotPCA(caf_dba_K36_norm, label=DBA_FACTOR)
dba.plotPCA(caf_dba_H4K20_norm, label=DBA_FACTOR)
dba.plotPCA(caf_dba_H3K4me2_norm, label=DBA_FACTOR)
dba.plotPCA(caf_dba_ATAC_norm, label=DBA_FACTOR)

#Venn Diagram of Gain vs loss compared to WT control
dba.plotVenn(caf_dba_K27_abc_norm, contrast = 4, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dba.plotVenn(caf_dba_K27_CS_norm, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dba.plotVenn(caf_dba_K36_norm, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dba.plotVenn(caf_dba_H4K20_norm, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dba.plotVenn(caf_dba_H3K4me2_norm, contrast = 3, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dba.plotVenn(caf_dba_ATAC_norm, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)

#MA Plots
dba.plotMA(caf_dba_K27_abc_norm, contrast = 4)
dba.plotMA(caf_dba_K27_CS_norm, contrast = 1)
dba.plotMA(caf_dba_K36_norm)
dba.plotMA(caf_dba_H4K20_norm)
dba.plotMA(caf_dba_H3K4me2_norm, contrast = 3)
dba.plotMA(caf_dba_ATAC_norm)

#Volcano Plots
dba.plotVolcano(caf_dba_K27_abc_norm, contrast = 1)
dba.plotVolcano(caf_dba_K27_CS_norm, contrast = 1)
dba.plotVolcano(caf_dba_K36_norm)
dba.plotVolcano(caf_dba_H4K20_norm)
dba.plotVolcano(caf_dba_H3K4me2_norm, contrast = 1)
dba.plotVolcano(caf_dba_ATAC_norm, contrast = 1)

#Box Plots
dba.plotBox(caf_dba_K27_abc_norm, contrast = 4)
dba.plotBox(caf_dba_K27_CS_norm, contrast = 1)
dba.plotBox(caf_dba_K36_norm)
dba.plotBox(caf_dba_H4K20_norm)
dba.plotBox(caf_dba_H3K4me2_norm, contrast = 3)
dba.plotBox(caf_dba_ATAC_norm, contrast = 1)

#Heatmaps
##For Binding Affinity
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(caf_dba_K27_abc_norm, contrast=1, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)
readscores <- dba.plotHeatmap(caf_dba_K27_CS_norm, contrast=1, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)
readscores <- dba.plotHeatmap(caf_dba_K36_norm, contrast=1, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)
readscores <- dba.plotHeatmap(caf_dba_H4K20_norm, contrast=1, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)
readscores <- dba.plotHeatmap(caf_dba_H3K4me2_norm, contrast=3, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)
readscores <- dba.plotHeatmap(caf_dba_ATAC_norm, contrast=1, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)

#Profiling
###system.file(’extra/plotProfileDemo.Rmd’,package=’DiffBind’)

#Basic
dba.plotProfile(caf_dba_K27_abc_norm, contrast = 1)
dba.plotProfile(caf_dba_K27_CS_norm, contrast = 1)
dba.plotProfile(caf_dba_K36_norm)
dba.plotProfile(caf_dba_H4K20_norm)
dba.plotProfile(caf_dba_H3K4me2_norm, contrast = 3)
dba.plotProfile(caf_dba_ATAC_norm, contrast = 1)

#split reps
mask.K27abc <- caf_dba$masks$abcam_H3K27me3
mask.K27CS <- caf_dba$masks$CS_H3K27me3
profiles <- dba.plotProfile(caf_dba, samples=list(H3K27me3_abcam= mask.K27abc, H3K27me3_CS= mask.K27CS), merge=NULL)
dba.plotProfile(profiles)

#Retrieve differentially bound sites for downstream analysis (will return a GRanges object)
## This is going to give me all sites that are different in each of the mutants vs. WT for each modification
caf_dba_K27_abc_norm.DB1 <- dba.report(caf_dba_K27_abc_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, bDB = TRUE)
caf_dba_K27_abc_norm.DB2 <- dba.report(caf_dba_K27_abc_norm, method=DBA_DESEQ2, contrast = 2, th=0.05,  bDB = TRUE)
caf_dba_K27_abc_norm.DB3 <- dba.report(caf_dba_K27_abc_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, bDB = TRUE)
caf_dba_K27_CS_norm.DB1 <- dba.report(caf_dba_K27_CS_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, bDB = TRUE)
caf_dba_K27_CS_norm.DB2 <- dba.report(caf_dba_K27_CS_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, bDB = TRUE)
caf_dba_K27_CS_norm.DB3 <- dba.report(caf_dba_K27_CS_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, bDB = TRUE)
caf_dba_K36_norm.DB1 <- dba.report(caf_dba_K36_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_K36_norm.DB2 <- dba.report(caf_dba_K36_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_K36_norm.DB3 <- dba.report(caf_dba_K36_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_K36_norm.DB4 <- dba.report(caf_dba_K36_norm, method=DBA_DESEQ2, contrast = 4, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_H4K20_norm.DB1 <- dba.report(caf_dba_H4K20_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, bDB = TRUE)
caf_dba_H4K20_norm.DB2 <- dba.report(caf_dba_H4K20_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, bDB = TRUE)
caf_dba_H4K20_norm.DB3 <- dba.report(caf_dba_H4K20_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, bDB = TRUE)
caf_dba_H3K4me2_norm.DB1 <- dba.report(caf_dba_H3K4me2_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_H3K4me2_norm.DB2 <- dba.report(caf_dba_H3K4me2_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_H3K4me2_norm.DB3 <- dba.report(caf_dba_H3K4me2_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_H3K4me2_norm.DB4 <- dba.report(caf_dba_H3K4me2_norm, method=DBA_DESEQ2, contrast = 4, th=0.05, fold = 0.8, bDB = TRUE)
caf_dba_ATAC_norm.DB1 <- dba.report(caf_dba_ATAC_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, bDB = TRUE)
caf_dba_ATAC_norm.DB2 <- dba.report(caf_dba_ATAC_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, bDB = TRUE)
caf_dba_ATAC_norm.DB3 <- dba.report(caf_dba_ATAC_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, bDB = TRUE)
caf_dba_ATAC_norm.DB4 <- dba.report(caf_dba_ATAC_norm, method=DBA_DESEQ2, contrast = 4, th=0.05, bDB = TRUE)
export.bed(caf_dba_ATAC_norm.DB3$peaks[[1]],"WT_v_cac-3_ATAC_DIFF.bed")

###saving reports
# Write to File
out <- as.data.frame(caf_dba_K27_abc_norm.DB1)
write.table(out, file="WT_v_cac-1_K27_CS_diff.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
bed <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(bed, file="WT_v_cac-3_ATAC_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


#########

#code to determine information that may be important/useful to your analysis

##Determine amount of peaks that are gained or lost in your samples, this will need to be tweaked depending on normalization method used
dbs <- dba.report(caf_dba, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
dbs$config$factor <- "normalize"
dbs$class[DBA_ID,] <- colnames(dbs$class)[1] <- "LIB_Full"
dbs$class[DBA_FACTOR,] <- DBA_NORM_LIB > dbs

##Compare FRiP scores to determine if ChIP efficency is cause of differential binding
dba.show(caf_dba,attributes=c(DBA_ID,DBA_FRIP))

##Plot a heatmap to compare the identified Differentialy bound peaks across normalization methods
deseq <- dba(dbs.all, mask=dbs.all$masks$DESeq2, minOverlap=1)
binding <- dba.peakset(deseq, bRetrieve=TRUE)
dba.plotHeatmap(deseq, maxSites=nrow(binding), bLog=FALSE, correlations=FALSE, minval=-5, maxval = 5, cexCol = 1.3, colScheme = hmap, main = "DESeq2 Differentially Bound Sites", ColAttributes = c(DBA_CONDITION, DBA_FACTOR), key.title = "LFC")

##Plot Peak Overlaps
olap.rate <- dba.overlap(caf_dba, mode=DBA_OLAP_RATE)
plot(olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')

##Genrate Consensus Peaksets for each modification analyzed
caf_consensus <- dba.peakset(caf_dba, consensus=-DBA_REPLICATE)
###combine consensus peaksets
caf_consensus <- dba(caf_consensus, mask=tamoxifen_consensus$masks$Consensus, minOverlap=1)

#########

#PLOTS

#Venn Diagram of Gain vs loss compared to WT control
caf_venn <- dba.plotVenn(caf_dba, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)

#PCA Plots
##All Peaks
dba.plotPCA(caf_dba,DBA_TREATMENT,label=DBA_FACTOR)

#Diff Bound Peaks
dba.plotPCA(caf_dba, contrast = 1,label=DBA_FACTOR)

#MA Plots
dba.plotMA(caf_dba)

#Volcano Plots
dab.plotVolcano(caf_dba)

#Box-Plots
pvals <- dba.plotBox(caf_dba)
plot(pvals)

#Heatmaps
##For Binding Affinity
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13) > readscores <- dba.plotHeatmap(caf_dba, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)

#Profiling
###system.file(’extra/plotProfileDemo.Rmd’,package=’DiffBind’)

#Basic
profiles <- dba.plotProfile(caf_dba)
dba.plotProfile(profiles)

#split reps
mask.K27abc <- caf_dba$masks$abcam_H3K27me3
mask.K27CS <- caf_dba$masks$CS_H3K27me3
profiles <- dba.plotProfile(caf_dba, samples=list(H3K27me3_abcam= mask.K27abc, H3K27me3_CS= mask.K27CS), merge=NULL)
dba.plotProfile(profiles)

#Area to write Multi-Factor designs (if motivated)

#Use this code to create a greylist of peaks enriched in inputs

#Use GreyListChIP

#caf_dba <- dba.blacklist(caf_dba)
#ChIP.greylist <- dba.blacklist(caf_dba, Retrieve=DBA_GREYLIST)



