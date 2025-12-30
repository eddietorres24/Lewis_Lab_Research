#commands to set up environment
#need to install DESeq2, ggplot2, and other common packages
#install scater package to calculate TPM from featurecounts
#uncomment the following three lines if you need to install the package

#local({r <- getOption("repos")
#r["CRAN"] <- "http://cran.r-project.org"
#options(repos=r)})

#if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")# BiocManager::install("scater")
# BiocManager::install("DESeq2")

library("DESeq2")
library("ggplot2")
library("gplots")
library("ggrepel")
library("dplyr")
library("pheatmap")
library("grid")
library("corrplot")
library(RColorBrewer)

##################################
#start anaysis here
##################################
#set working directory to folder containing featureCounts output
setwd("C:\\Users\\eddie\\Desktop\\Research\\GitHub")

#Bring in table with unnormalized transcription counts; check.names is important if you have dashes in the gene names 
#row.names=1 command sets geneIDs as the row name
countdataInteractors <- read.table("./readcounts_FINAL.txt",skip=1, header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
#countdataInteractors_ash3 <- read.table("./readcounts_ash3.txt",skip=1, header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
countdataInteractors_final = cbind(countdataInteractors, countdataInteractors_ash3[,6:7])
countdataInteractors <- data.matrix(countdataInteractors_final)

#################################################
#Section 1: calculate tpm using the scater package
#################################################

#convert count file to matrix. Note:: I initially tried to use as.matrix, 
#but this DID NOT WORK. once I used v <- data.matrix(countfile) then it all worked fine#

#create matrix containing count data only; i.e. eliminate contig, start, stop, length fields
Interactors_countsOnly <- data.matrix(countdataInteractors_final[ ,6:ncol(countdataInteractors_final)])

##Run calculateRPM function - must load "scater" library
library("scater")
Interactors_tpm <- calculateTPM(Interactors_countsOnly, lengths = countdataInteractors_final[,5])

##convert to matrix again
allDataTPM  <- data.matrix(Interactors_tpm)

#rename columns with human readable sample names

###THIS NEXT SECTION REQUIRED MANUAL INPUT INTO A SPREADSHEET OR TEXT EDITOR

#step 1 create a varialbe with the current, long column names
samplesname <- colnames(allDataTPM)

#write the long sample names to a spreadsheet file so these can be manipulated in excel or text editor
write.table(samplesname, file="samplenames.txt", sep="\t")

#################

#create a list of column names that you want to plot; !!!!You must include the column name for gene ID

#move the subset of genes you want to plot into a new matrix

sampleNames <- colnames(allDataTPM)
write.table(sampleNames, file="stuff.txt")
Ordered_KO_data <- cbind(allDataTPM[,22:25],allDataTPM[,26:28],allDataTPM[,1:3],allDataTPM[,4:6],allDataTPM[,7:9],allDataTPM[,10:12],allDataTPM[,13:15],allDataTPM[,16:18],allDataTPM[,19:21],allDataTPM[,29:30])
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,22:25], na.rm = TRUE),
                                 rowMeans(allDataTPM[,26:28], na.rm = TRUE),
                                 rowMeans(allDataTPM[,29:30], na.rm = TRUE),
                                 rowMeans(allDataTPM[,1:3], na.rm = TRUE),
                                 rowMeans(allDataTPM[,4:6], na.rm = TRUE),
                                 rowMeans(allDataTPM[,7:9], na.rm = TRUE),
                                 rowMeans(allDataTPM[,10:12], na.rm = TRUE),
                                 rowMeans(allDataTPM[,13:15], na.rm = TRUE),
                                 rowMeans(allDataTPM[,16:18], na.rm = TRUE),
                                 rowMeans(allDataTPM[,19:21], na.rm = TRUE))
averageRowIDs=c("WT","set-7","set-3","cac-1","cac-2","cac-3","naf-1","asf-1","naf-2","ATRX")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

####################################################################

#SECTION 3
#subset data to sample only PRC2-target genes

#read in geneIDs of PRC2 target genes
Prc2targets <- read.table("heatmapPRC2genesEdit.txt", header=FALSE,stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

###SUBSET DATA FOR ALL KO SAMPLES
Prc2targetTPM <- subset(allDataTPM, !(rownames(allDataTPM)%in%Prc2targets[,1]))
Prc2targetTPM <- subset(Prc2targetTPM, rownames(Prc2targetTPM)%in%rownames(cac1DE) | rownames(Prc2targetTPM)%in%rownames(cac2DE) | rownames(Prc2targetTPM)%in%rownames(cac3DE))

###SUBSET DATA FOR AVERAGED KO SAMPLES
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%Prc2targets[,1])

###SUBSET TO GET RID OF PRC2 GENES, BUT KEEP OTHER DE GENES (cac-1/2/3)
AVERAGE_NonPRC2targetTPM <- subset(Averaged_Orderd_KO_data, !(rownames(Averaged_Orderd_KO_data)%in%Prc2targets[,1]))
AVERAGE_NonPRC2targetTPM <- subset(AVERAGE_NonPRC2targetTPM, rownames(AVERAGE_NonPRC2targetTPM)%in%rownames(cac1DE) | rownames(AVERAGE_NonPRC2targetTPM)%in%rownames(cac2DE) | rownames(AVERAGE_NonPRC2targetTPM)%in%rownames(cac3DE))

###Subset data to filter out non-PRC2 target regions that got through (cutting out any genes over 2 tpm in WT)
AVERAGE_Prc2targetTPM <- subset(AVERAGE_Prc2targetTPM, (AVERAGE_Prc2targetTPM[,1] < 2))

###Add sudocount and log transform (if necessary)
AVERAGE_Prc2targetTPM <- AVERAGE_Prc2targetTPM + 0.1
AVERAGE_Prc2targetTPM <- log2(AVERAGE_Prc2targetTPM)

#Non-PRC2
AVERAGE_NonPRC2targetTPM <- AVERAGE_NonPRC2targetTPM + 0.1
AVERAGE_NonPRC2targetTPM <- log2(AVERAGE_NonPRC2targetTPM)

#melt data to get it into a format ggplots can use (load library "reshape2")
library(reshape2)

meltedPRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

#################################################################################

##PLOT AVERAGE DATA FOR FIGURE 1A

##################################################################
#melt the data and label value "Count" and the sample IDs GeneID and Sample
#meltedData <- melt(tpm, value.name = "Count",
#varnames=c('GeneID', 'Sample'))

meltedAveragePRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

altorder = rev(c( "WT","set-7","set-3","cac-1","cac-2","cac-3","naf-1","naf-2","asf-1","ATRX"))
meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample, altorder)

# Plot box & whisker chart
library(ggplot2)
xlabels = averageRowIDs
colors = rev(c( "#4575b4","#fee090","#fee090","#fee090", "#fee090", "#fee090", "#4575b4", "#4575b4", "#4575b4", "#4575b4"))

box<-ggplot(meltedAveragePRC2targetData, aes(x=Sample, y=Count)) +
  labs(y="Expression Level (Transcripts per Million)", x="Strain") +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill=colors, size=0.1, coef=1.5, lwd=0.25) +
  coord_flip(ylim=c(0,15))+
  theme_get() + 
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=3, colour = "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  theme(text = element_text(size = 10), axis.title.x=element_text(size=10, margin=margin(t = 10, r = 0, b = 0, l = 0),),axis.title.y=element_text( size=10)) +
  theme(axis.text.x=element_text(angle = 270, vjust = 0.5, hjust = 1, size=10 )) +
  theme(axis.text.y=element_text(vjust = 0.5, hjust = 1, size=10 )) 
  guides(fill=F) 

#device.on()
box

ggsave(filename = "histone_chaperone_boxplot_PAPER.pdf", plot = box, dpi=600, height= 3, width=4)


########################3
############################################################################

library(pheatmap)
library(RColorBrewer)

breaks1=seq(-4, 5, by=.09) #This is to set a custom heatmaps scale. Not used here.

#filter all data with no changes >>>>> if the sum of all rows divided by the # of rows is > 0.

##Experiment with different heatmap strategies
##strategy 1 - Plot the average TMP of all genes that change expression; row normalization; clustredRows

##Keep Only Genes that are expressed in at least one sample
GenesWithChanges <- subset(AVERAGE_NonPRC2targetTPM, (rowSums(allDataTPM) > 0))
GenesWithChanges_95 <- subset(AVERAGE_Prc2targetTPM, (rowSums(Prc2targetTPM) > 0))

#plot in the desired column order; did this by subsetting the dataset based on sample list 'altorder' above

heatmap<- pheatmap(GenesWithChanges_no_PRC2[,rev(altorder)], color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                   cellwidth = NA, cellheight = NA, cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                   legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot<-heatmap[[4]]

ggsave(filename = "./histone_chaperone_heatmap_log2_prc2_gene_no_scale.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)
#dev.off()