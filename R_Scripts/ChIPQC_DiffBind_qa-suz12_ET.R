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
qasamples <- read.csv("./DiffBind_qa-suz12_fig4_K27_ET.csv")


## Make ChIPQC Object
register(SerialParam())
chipObj <- ChIPQC(qasamples) 

## Make QC Report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_qa-suz12_fig4_H3K27me3", reportFolder="ChIPQCreport_qa-suz12")

### DiffBind (scratch)

#Read in sample sheet containing CAF-1 mutant data *If running ChIP QC, you can use "samples" variable
qa_dba <- dba(sampleSheet = qasamples)


##Plot correlation plot between all samples included on data sheet, save in a file
CorrPlot_qa <- plot(qa_dba)

#Using RPKM
CorrPlot_RPKM <- dba.plotHeatmap(qa_dba, score=DBA_SCORE_RPKM_FOLD)

#Count Reads
qa_dba <- dba.count(qa_dba)

##Plot new correlation plot based on count data, save in a file
CorrPlot_count_qa <- plot(qa_dba)

#Normalization
qa_dba_norm <- dba.normalize(qa_dba, normalize="lib")

#Model design & Contrast (what comparisons do you want to make?)
qa_dba_norm <- dba.contrast(qa_dba_norm, categories=DBA_FACTOR, minMembers = 2)

#Blacklist
#qa_dba <- dba.blacklist(qa_dba, blacklist=FALSE, greylist=FALSE)

#Differential Analysis
qa_dba_norm <- dba.analyze(qa_dba_norm)

#Plot PCA for All mods (normalized peak calls)
dba.plotPCA(qa_dba_norm, label=DBA_FACTOR)

#Venn Diagram of Gain vs loss compared to WT control
dba.plotVenn(qa_dba_norm, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)

#MA Plots
dba.plotMA(qa_dba_norm, contrast = 1)

#Volcano Plots
dba.plotVolcano(qa_dba_norm, contrast = 4)

#Box Plots
dba.plotBox(qa_dba_norm, contrast = 1)

#Heatmaps
##For Binding Affinity
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(qa_dba_norm, contrast=4, correlations=FALSE, scale="row", colScheme = hmap, cexCol = 0.9)

#Profiling
###system.file(’extra/plotProfileDemo.Rmd’,package=’DiffBind’)

#Basic
dba.plotProfile(qa_dba_norm, contrast = 1)

#split reps
#mask.K27abc <- caf_dba$masks$abcam_H3K27me3
#mask.K27CS <- caf_dba$masks$CS_H3K27me3
#profiles <- dba.plotProfile(caf_dba, samples=list(H3K27me3_abcam= mask.K27abc, H3K27me3_CS= mask.K27CS), merge=NULL)
#dba.plotProfile(profiles)

#Retrieve differentially bound sites for downstream analysis (will return a GRanges object)
## This is going to give me all sites that are different in each of the mutants vs. WT for each modification
qa_dba_norm.DB1 <- dba.report(qa_dba_norm, method=DBA_DESEQ2, contrast = 1, th=0.05, bDB = TRUE)
qa_dba_norm.DB2 <- dba.report(qa_dba_norm, method=DBA_DESEQ2, contrast = 2, th=0.05, bDB = TRUE)
qa_dba_norm.DB3 <- dba.report(qa_dba_norm, method=DBA_DESEQ2, contrast = 3, th=0.05, bDB = TRUE)
qa_dba_norm.DB4 <- dba.report(qa_dba_norm, method=DBA_DESEQ2, contrast = 4, th=0.05, bDB = TRUE)
qa_dba_norm.DB5 <- dba.report(qa_dba_norm, method=DBA_DESEQ2, contrast = 5, th=0.05, bDB = TRUE)

###saving reports
export.bed(qa_dba_norm.DB1$peaks[[1]],"WT_v_0hr_DIFF.bed")

# Write to File
out <- as.data.frame(caf_dba_K27_abc_norm.DB1)
write.table(out, file="WT_v_cac-1_K27_CS_diff.txt", sep="\t", quote=F, row.names=F)

#old way
# Create bed files for each keeping only significant peaks (p < 0.05)
#bed <- out %>% 
#  filter(FDR < 0.05 & Fold > 0) %>% 
#  select(seqnames, start, end)

# Write to file
#write.table(bed, file="WT_v_cac-3_ATAC_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


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



