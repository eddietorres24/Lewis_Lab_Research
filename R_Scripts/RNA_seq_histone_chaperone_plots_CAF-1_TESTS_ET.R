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

knitr::opts_chunk$set(echo = TRUE)

# Set working environment

workingdir="C:/Users/eddie/Research/GitHub"

#set working directory to the correct location for working machine
knitr::opts_knit$set(root.dir = "workingdir")

##################################
#start anaysis here
##################################

#Bring in table with un-normalized transcription counts; check.names is important if you have dashes in the gene names 
#row.names=1 command sets geneIDs as the row name
countdataInteractors <- read.table("./text_files/readcounts_FINAL.txt",skip=1, header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")

#################################################
#Section 1: calculate tpm using the scater package
#################################################

#convert count file to matrix. Note:: I initially tried to use as.matrix, 
#but this DID NOT WORK. once I used v <- data.matrix(countfile) then it all worked fine#

#create matrix containing count data only; i.e. eliminate contig, start, stop, length fields
Interactors_countsOnly <- data.matrix(countdataInteractors[ ,6:ncol(countdataInteractors)])

##Run calculateTPM function - must load "scater" library
library("scater")
Interactors_tpm <- calculateTPM(Interactors_countsOnly, lengths = countdataInteractors[,5])

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
write.table(sampleNames, file="names.txt")
Ordered_KO_data <- cbind(allDataTPM[,22:25],allDataTPM[,26:28],allDataTPM[,29:31],allDataTPM[,4:6],allDataTPM[,32:34],allDataTPM[,7:9],allDataTPM[,10:12],allDataTPM[,16:18],allDataTPM[,13:15],allDataTPM[,19:21])
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,22:25], na.rm = TRUE),
                                 rowMeans(allDataTPM[,26:28], na.rm = TRUE),
                                 rowMeans(allDataTPM[,29:31], na.rm = TRUE),
                                 rowMeans(allDataTPM[,4:6], na.rm = TRUE),
                                 rowMeans(allDataTPM[,32:34], na.rm = TRUE),
                                 rowMeans(allDataTPM[,7:9], na.rm = TRUE),
                                 rowMeans(allDataTPM[,10:12], na.rm = TRUE),
                                 rowMeans(allDataTPM[,16:18], na.rm = TRUE),
                                 rowMeans(allDataTPM[,13:15], na.rm = TRUE),
                                 rowMeans(allDataTPM[,19:21], na.rm = TRUE))
averageRowIDs=c("WT","set-7","cac-1","cac-2","cac-1_cac-2","cac-3","naf-1","naf-2","asf-1","ATRX")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

####################################################################

#SECTION 3
#subset data to sample only PRC2-target genes (or other gene sets of interest)

#read in geneIDs of PRC2 target genes
#reading in csvs w/ upregulated genes in CAF-1 mutants
Prc2targets <- read.table("./bed_files/K27_genes_gff_sorted_names.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 
cac1up <- read.csv("./cac_DEseq/cac1_UP.csv", stringsAsFactors=FALSE, check.names=FALSE)
cac2up <- read.csv("./cac_DEseq/cac2_UP.csv", stringsAsFactors=FALSE, check.names=FALSE)

###SUBSET DATA FOR ALL KO SAMPLES
Prc2targetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%Prc2targets[,1])

#CAF-1 Upregulated mutants
CAFUpTPM <- subset(allDataTPM, rownames(allDataTPM) %in% cac1up$NCU | rownames(allDataTPM) %in% cac2up$NCU)

###SUBSET DATA FOR AVERAGED KO SAMPLES
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%Prc2targets[,1])
AVERAGE_CAFTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data) %in% cac1up$NCU | rownames(allDataTPM) %in% cac2up$NCU)

###Subset data to filter out non-PRC2 target regions that got through (cutting out any genes over 2 tpm in WT)
AVERAGE_Prc2targetTPM <- subset(AVERAGE_Prc2targetTPM, (AVERAGE_Prc2targetTPM[,1] < 12.5))
AVERAGE_AlldataTPM <- subset(Averaged_Orderd_KO_data, (Averaged_Orderd_KO_data[,1] > -0.1))

###resubet PRC2 targets after filtering
Prc2targetTPM <- subset(Prc2targetTPM, rownames(Prc2targetTPM)%in%rownames(AVERAGE_Prc2targetTPM))

###Add sudocount and log transform (if necessary)
AVERAGE_Prc2targetTPM <- AVERAGE_Prc2targetTPM + 1.1
AVERAGE_Prc2targetTPM <- log2(AVERAGE_Prc2targetTPM)

AVERAGE_CAFTPM <- AVERAGE_CAFTPM + 1
AVERAGE_CAFTPM <- log2(AVERAGE_CAFTPM)

AVERAGE_AlldataTPM <- AVERAGE_AlldataTPM + 1
AVERAGE_AlldataTPM <- log2(AVERAGE_AlldataTPM)
#melt data to get it into a format ggplots can use (load library "reshape2")
library(reshape2)

meltedPRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

#Alldata
meltedAllData <- melt(AVERAGE_AlldataTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))


#################################################################################

##PLOT AVERAGE DATA FOR FIGURE 1A

##################################################################
#melt the data and label value "Count" and the sample IDs GeneID and Sample
#meltedData <- melt(tpm, value.name = "Count",
#varnames=c('GeneID', 'Sample'))

meltedAveragePRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

meltedNotInCac <- melt(K27_lost, value.name = 'Count',
                                    varnames=c('GeneID', 'Sample'))

meltedInCac <- melt(K27_not_lost, value.name = 'Count',
                                    varnames=c('GeneID', 'Sample'))

meltedAverageAllData <- melt(AVERAGE_AlldataTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

altorder = rev(c("WT","set7","cac1","cac2","cac1cac2","cac3","naf1","naf2","asf1","ATRX"))
altorder1 = rev(c("WT","1"))

meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample)

meltedAverageAllData$Sample <- factor(meltedAverageAllData$Sample)

# Plot box & whisker chart
library(ggplot2)
xlabels = averageRowIDs
colors = c( "#4575b4","#fee090","#fee090", "#fee090", "#fee090", "#fee090","#4575b4", "#4575b4", "#4575b4", "#4575b4")

box<-ggplot(meltedAveragePRC2targetData, aes(x=Sample, y=Count)) +
  labs(y="Expression Level (Transcripts per Million)", x="Strain") +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(notch = TRUE, fill=colors, size=0.1, coef=1.5, lwd=0.25) +
  ylim(0, 11.5) +
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

ggsave(filename = "histone_chaperone_boxplot_PAPER_FINAL.pdf", plot = box, dpi=600, height= 3, width=4)

#######################

#t-test for significance in difference of mean expression values of PRC2 targets b/w WT & mutants

set7_t <- t.test(AVERAGE_Prc2targetTPM[,2], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
cac1_t <- t.test(AVERAGE_Prc2targetTPM[,3], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
cac2_t <- t.test(AVERAGE_Prc2targetTPM[,4], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
cac1_2_t <- t.test(AVERAGE_Prc2targetTPM[,5], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
cac3_t <- t.test(AVERAGE_Prc2targetTPM[,6], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
naf1_t <- t.test(AVERAGE_Prc2targetTPM[,7], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
naf2_t <- t.test(AVERAGE_Prc2targetTPM[,8], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
asf1_t <- t.test(AVERAGE_Prc2targetTPM[,9], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)
atrx_t <- t.test(AVERAGE_Prc2targetTPM[,10], AVERAGE_Prc2targetTPM[,1], data = AVERAGE_Prc2targetTPM)

#I want to make a violin plot instead of boxplot. Going to adapt Abby's code
library(grDevices)

total2 <- meltedInCac %>%
  group_by(Sample)
total_dist = meltedInCac %>%
  group_by(Sample) %>% summarise(num=n())

#Setting distance between violins in plot
dodge <- position_dodge(width = 1)
###  plotting violin plot
### this chunk is just for putting samples in the order that I want ###
violin <- total2 %>%
  left_join(total_dist) %>%
  arrange(factor(Sample, levels = c("WT", "set-7", "cac-1", "cac-1_new", "cac-2", "cac-1_cac-2", "cac-1_suz12", "cac-3", "naf-1", "naf-2", "asf-1", "ATRX"))) %>%
  mutate(Sample = factor(Sample)) %>%
  ### everything below is the actual violin plot ###
  ggplot(aes(x=Sample, y=Count)) + 
  geom_violin(position = dodge, scale="width", trim=FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
  scale_fill_manual("",values = c("orchid1", "springgreen3")) +
  labs(x = "Strain",y = expression("Expression Level (log "[2]~"(TPM+1)")) + 
  theme_classic(base_size = 20)

print(violin)
ggsave("./K27_genes_In_Cac.pdf", plot=violin, width = 10, height = 8, unit="in",  dpi=400)

#stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, colour = "red") +

########################3
############################################################################

library(pheatmap)
library(RColorBrewer)

breaks1=seq(-4, 5, by=.09) #This is to set a custom heatmaps scale. Not used here.

#filter all data with no changes >>>>> if the sum of all rows divided by the # of rows is > 0.

##Experiment with different heatmap strategies
##strategy 1 - Plot the average TMP of all genes that change expression; row normalization; clustredRows

##Keep Only Genes that are expressed in at least one sample
GenesWithChanges <- subset(AVERAGE_Prc2targetTPM, (rowSums(AVERAGE_Prc2targetTPM) > 0))
#CAFUpGenes <- subset(AVERAGE_CAFTPM, (rowSums(CAFUpTPM) > 0))
#GenesWithChanges_95 <- subset(AVERAGE_Prc2targetTPM, (rowSums(Prc2targetTPM) > 0))
K27_lost <- subset(AVERAGE_Prc2targetTPM, rownames(AVERAGE_Prc2targetTPM) %in% K27_genes_notincac$V10)
K27_lost <- subset(K27_lost, rowSums(K27_lost) > 0)
K27_not_lost <- subset(AVERAGE_Prc2targetTPM, ! rownames(AVERAGE_Prc2targetTPM) %in% K27_genes_notincac$V10)
K27_not_lost <- subset(K27_not_lost, rowSums(K27_not_lost) > 0)

#plot in the desired column order; did this by subsetting the dataset based on sample list 'altorder' above

heatmap <- pheatmap(GenesWithChanges, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                   cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                   legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

heatmap <- pheatmap(combine, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                    cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = F, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                    legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot <- heatmap[[4]]

ggsave(filename = "./CAF-1_K27_lost.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)
#dev.off()