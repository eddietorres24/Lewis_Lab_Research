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
library(reshape2)
library(grDevices)
library("scater")

knitr::opts_chunk$set(echo = TRUE)

# Set working environment

workingdir="C:/Users/eddie/Research/GitHub"

#set working directory to the correct location for working machine
knitr::opts_knit$set(root.dir = "workingdir")

##################################
#start analysis here, qa-suz12
##################################


#Bring in table with unnormalized transcription counts; check.names is important if you have dashes in the gene names 
#row.names=1 command sets geneIDs as the row name
countdataInteractors <- read.table("./text_files/readcounts_qa_paper.txt",skip=1, header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
#countdataInteractors_final = cbind(countdataInteractors, countdataInteractors_ash3[,6:7])
#countdataInteractors <- data.matrix(countdataInteractors_final)

#################################################
#Section 1: calculate tpm using the scater package
#################################################

#convert count file to matrix. Note:: I initially tried to use as.matrix, 
#but this DID NOT WORK. once I used v <- data.matrix(countfile) then it all worked fine#

#create matrix containing count data only; i.e. eliminate contig, start, stop, length fields
Interactors_countsOnly <- data.matrix(countdataInteractors[ ,6:ncol(countdataInteractors)])

##Run calculateTPM function - must load "scater" library
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
Ordered_KO_data <- cbind(allDataTPM[,4:6],allDataTPM[,7:9],allDataTPM[,1:3],allDataTPM[,10:12],allDataTPM[,13:15])
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,4:6], na.rm = TRUE),
                                 rowMeans(allDataTPM[,7:9], na.rm = TRUE),
                                 rowMeans(allDataTPM[,1:3], na.rm = TRUE),
                                 rowMeans(allDataTPM[,10:12], na.rm = TRUE),
                                 rowMeans(allDataTPM[,13:15], na.rm = TRUE))
averageRowIDs=c("WT_0hr","WT_24hr","suz-12","qa-suz12_0hr", "qa-suz12_24hr")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs


####################################################################

#SECTION 3
#subset data to sample only PRC2-target genes

#read in geneIDs of PRC2 target genes
Prc2targets <- read.table("./bed_files/H3K27me3_methylated_genes_FINAL.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

#reading in csvs w/ upregulated genes in CAF-1 mutants

###SUBSET DATA FOR ALL KO SAMPLES
Prc2targetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%Prc2targets[,10])

###SUBSET DATA FOR AVERAGED KO SAMPLES
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%Prc2targets[,10])

###Subset data to filter out non-PRC2 target regions that got through (cutting out any genes over 2 tpm in WT)
AVERAGE_Prc2targetTPM <- subset(AVERAGE_Prc2targetTPM, (AVERAGE_Prc2targetTPM[,1] < 5))
AVERAGE_AlldataTPM <- subset(Averaged_Orderd_KO_data, (Averaged_Orderd_KO_data[,1] > -0.1))

###resubet PRC2 targets after filtering
Prc2targetTPM <- subset(Prc2targetTPM, rownames(Prc2targetTPM)%in%rownames(AVERAGE_Prc2targetTPM))

###Add sudocount and log transform (if necessary)
AVERAGE_Prc2targetTPM <- AVERAGE_Prc2targetTPM + 1
AVERAGE_Prc2targetTPM <- log2(AVERAGE_Prc2targetTPM)

AVERAGE_AlldataTPM <- AVERAGE_AlldataTPM + 1
AVERAGE_AlldataTPM <- log2(AVERAGE_AlldataTPM)
#melt data to get it into a format ggplots can use (load library "reshape2")


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

meltedAverageAllData <- melt(AVERAGE_AlldataTPM, value.name = 'Count',
                      varnames=c('GeneID', 'Sample'))

altorder = rev(c("WT_0hr","WT_24hr","suz-12","qa-suz12_0hr", "qa-suz12_24hr"))

meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample)

meltedAverageAllData$Sample <- factor(meltedAverageAllData$Sample)

# Plot box & whisker chart
xlabels = averageRowIDs
colors = c( "#4575b4","#4575b4", "#fee090", "#fee090", "orange")

box<-ggplot(meltedAveragePRC2targetData, aes(x=Sample, y=Count)) +
  labs(y="Expression Level (Transcripts per Million)", x="Strain") +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(notch = TRUE, fill=colors, size=0.1, coef=1.5, lwd=0.25) +
  ylim(0, 10) +
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

ggsave(filename = "histone_chaperone_boxplot_PAPER_DELETE.pdf", plot = box, dpi=600, height= 3, width=4)

#######################

#t-test for significance in difference of mean expression values of PRC2 targets b/w WT & mutants

# set7_t <- t.test(PRC2meancalc$set.7, mu = mean(PRC2meancalc$WT))
# cac1_t <- t.test(PRC2meancalc$cac.1, mu = mean(PRC2meancalc$WT))
# cac2_t <- t.test(PRC2meancalc$cac.2, mu = mean(PRC2meancalc$WT))
# cac3_t <- t.test(PRC2meancalc$cac.3, mu = mean(PRC2meancalc$WT))
# naf1_t <- t.test(PRC2meancalc$naf.1, mu = mean(PRC2meancalc$WT))
# naf2_t <- t.test(PRC2meancalc$naf.2, mu = mean(PRC2meancalc$WT))
# asf1_t <- t.test(PRC2meancalc$asf.1, mu = mean(PRC2meancalc$WT))
# atrx_t <- t.test(PRC2meancalc$ATRX, mu = mean(PRC2meancalc$WT))

#I want to make a violin plot instead of boxplot. Going to adapt Abby's code


total2 <- meltedAveragePRC2targetData %>%
  group_by(Sample)
total_dist = meltedAveragePRC2targetData %>%
  group_by(Sample) %>% summarise(num=n())

#Setting distance between violins in plot
dodge <- position_dodge(width = 1)
###  plotting violin plot
### this chunk is just for putting samples in the order that I want ###
violin <- total2 %>%
  left_join(total_dist) %>%
  arrange(factor(Sample, levels = c("WT_0hr","WT_24hr","suz-12","qa-suz12_0hr", "qa-suz12_24hr"))) %>%
  mutate(Sample = factor(Sample)) %>%
  ### everything below is the actual violin plot ###
  ggplot(aes(x=Sample, y=Count)) + 
  geom_violin(position = dodge, scale="width", trim=FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
  scale_fill_manual("",values = c("orchid1", "springgreen3")) +
  labs(x = "Strain",y = expression("Expression Level (log "[2]~"(TPM+1)")) + 
  theme_classic(base_size = 20)

print(violin)
ggsave("./qa-suz12_PRC2_Gene_expression.pdf", plot=violin, width = 10, height = 8, unit="in",  dpi=400)

#stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, colour = "red") +

########################3
############################################################################

breaks1=seq(-4, 5, by=.09) #This is to set a custom heatmaps scale. Not used here.

#filter all data with no changes >>>>> if the sum of all rows divided by the # of rows is > 0.

##Experiment with different heatmap strategies
##strategy 1 - Plot the average TMP of all genes that change expression; row normalization; clustredRows

##Keep Only Genes that are expressed in at least one sample
GenesWithChanges <- subset(AVERAGE_Prc2targetTPM, (rowSums(Prc2targetTPM) > 0) & (rowSums(AVERAGE_Prc2targetTPM) > 0))
#CAFUpGenes <- subset(AVERAGE_CAFTPM, (rowSums(CAFUpTPM) > 0))
GenesWithChanges_95 <- subset(AVERAGE_Prc2targetTPM, (rowSums(Prc2targetTPM) > 0))

#plot in the desired column order; did this by subsetting the dataset based on sample list 'altorder' above

heatmap1 <- pheatmap(GenesWithChanges[,rev(altorder)], color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                   cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                   legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot <- heatmap1[[4]]

ggsave(filename = "./qa-suz12_log2_prc2_genes_heatmap_rowscale.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)
#dev.off()
