### Script for ChIPQC


## Load Libraries
library(BiocManager)
library(ChIPQC)
library(ggplot2)

## Load Sample Data
setwd("C:\\Users\\eddie\\Research\\GitHub")
samples <- read.csv("./DiffBind_CAF-1_ET.csv")

## Make ChIPQC Object
register(SerialParam())
chipObj <- ChIPQC(samples) 

## Make QC Report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_CAF-1_mods", reportFolder="ChIPQCreport")


### DiffBind (scratch)

#Read in sample sheet containing CAF-1 mutant data *If running ChIP QC, you can use "samples" variable
caf_dba <- dba(sampleSheet = samples)

##Plot correlation plot between all samples included on data sheet, save in a file
CorrPlot <- plot(caf_dba)
ggsave(filename = "CAF-1_DBA_Corr_Plot.png")
#Using RPKM
CorrPlot_RPKM <- dba.plotHeatmap(caf_dba, score=DBA_SCORE_RPKM_FOLD)

#Count Reads
caf_dba <- dba.count(caf_dba)

##Plot new correlation plot based on count data, save in a file
CorrPlot_count <- plot(caf_dba)
ggsave(filename = "CAF-1_DBA_Corr_Plot_counts.pdf", plot = last_plot(), dpi=600)

#Normalization
caf_dba <- dba.normalize(caf_dba, method = DBA_ALL_METHODS, background = TRUE)

#Model design & Contrast (what comparisons do you want to make?)
caf_dba <- dba.contrast(caf_dba, design = "~Factor + Treatment")

#Differential Analysis
caf_dba <- dba.analyze(caf_dba)

##Plot correlation plot based only on differentially bound sites, save in file
CorrPlot_diff <- plot(caf_dba, contrast = 1)
ggsave(filename = "CAF-1_DBA_Corr_Plot_diff.pdf", plot = CorrPlot_diff, dpi=600)

#Retrieve differentially bound sites for downstream analysis (will return a GRanges object)
caf_dba.DB <- dba.report(caf_dba)


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



### DiffBind Script

## Load Libraries
library(DiffBind)
library(tidyverse)

## Run DBA
dbObj <- dba(sampleSheet=samples)

## Plot PCA & Correlation Plot
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)

## Contrast Samples
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)

## Differential Enrichment Analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)

# Plot
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)

## Visualize Results

# Peaks recognized by edgeR vs. DEseq
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

# MA Plots
dba.plotMA(dbObj, contrast=1, method=DBA_DESEQ2)
dba.plotMA(dbObj, contrast=1, bXY=TRUE)
pvals <- dba.plotBox(dbObj, contrast=1)

## Extract Results
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 3, th=1)

# Write to File
out <- as.data.frame(res_deseq)
write.table(out, file="cac3_H3K27me3_DiffBind.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)

WT_cac3_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(WT_cac3_enrich, file="WT_cac3_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

cac3_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(cac3_enrich, file="cac3_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

