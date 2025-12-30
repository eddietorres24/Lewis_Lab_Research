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
#library("ggpubr")
library(RColorBrewer)


##################################
#start anaysis here
##################################
#set working directory to folder containing featureCounts output
setwd("C:\\Users\\eddie\\Desktop\\Research\\GitHub")

#Bring in table with unnormalized transcription counts; check.names is important if you have dashes in the gene names 
#row.names=1 command sets geneIDs as the row name
countdataInteractors <- read.table("./readcounts_no_isw_rtt109.txt",skip=1, header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")

countdataInteractors <- data.matrix(countdataInteractors)


#################################################
#Section 1: calculate tpm using the scater package
#################################################

#convert count file to matrix. Note:: I initially tried to use as.matrix, 
#but this DID NOT WORK. once I used v <- data.matrix(countfile) then it all worked fine#

#create matrix containing count data only; i.e. eliminate contig, start, stop, length fields
Interactors_countsOnly<- data.matrix(countdataInteractors[ ,6:ncol(countdataInteractors)])

##Run calculateRPM function - must load "scater" library

library("scater")
Interactors_tpm <- calculateTPM(Interactors_countsOnly, lengths = countdataInteractors[,5])

##convert to matrix again
allDataTPM  <- data.matrix(Interactors_tpm)

####################
#rename columns with human readable sample names

###THIS NEXT SECTION REQUIRED MANUAL INPUT INTO A SPREADSHEET OR TEXT EDITOR

#step 1 create a varialbe with the current, long column names
samplesname <- colnames(allDataTPM)

#write the long sample names to a spreadsheet file so these can be manipulated in excel or text editor
write.table(samplesname, file="samplenames.txt", sep="\t")

#################
#Section 2: Skip this section for now_this is for percentile normalizing for heat mapping, 

#normalizaton methods; do this after it is working

#calculate the 95th percentile value
####max = quantile(tpm,0.95)

#replace outliers with 95th percentile value
####tpm[tpm>max] = max

#remove rows with zeros for every gene; these cause white lines to be produced in the heatmap
###data_mat <- data_mat[apply(data_mat[,-1], 1, function(x) !all(x==0)),]
###data_mat<-as.matrix(data_mat)

#add row names
####rownames(data_mat) <- dataRnames  

#create a list of column names that you want to plot; !!!!You must include the column name for gene ID

#move the subset of genes you want to plot into a new matrix

sampleNames <- colnames(allDataTPM)
write.table(sampleNames, file="stuff.txt")
Ordered_KO_data <- cbind(allDataTPM[,22:25],allDataTPM[,26:28],allDataTPM[,1:3],allDataTPM[,4:6],allDataTPM[,7:9],allDataTPM[,10:12],allDataTPM[,13:15],allDataTPM[,16:18],allDataTPM[,19:21])
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,22:25], na.rm = TRUE),
                                 rowMeans(allDataTPM[,26:28], na.rm = TRUE),
                                 rowMeans(allDataTPM[,1:3], na.rm = TRUE),
                                 rowMeans(allDataTPM[,4:6], na.rm = TRUE),
                                 rowMeans(allDataTPM[,7:9], na.rm = TRUE),
                                 rowMeans(allDataTPM[,10:12], na.rm = TRUE),
                                 rowMeans(allDataTPM[,13:15], na.rm = TRUE),
                                 rowMeans(allDataTPM[,16:18], na.rm = TRUE),
                                 rowMeans(allDataTPM[,19:21], na.rm = TRUE))
averageRowIDs=c("WT","set-7","cac-1","cac-2","cac-3","naf-1","asf-1","naf-2","ATRX")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

###calculate summary stats for all gene data
TPMmean <- colMeans(Averaged_Orderd_KO_data)
TPMmedian <- colMedians(Averaged_Orderd_KO_data)
TPMstdev <- apply(Averaged_Orderd_KO_data,2,sd)
w<- cbind(TPMmean,TPMmedian,TPMstdev)

#write.table(w, file="~/Dropbox/DropBOX Documents/Zack Papers/2018_RemodellerPaper/IswiPaperFilesFromMasayuki/SummaryTablesFromGenomics/KO_AverageTPM_stats.txt", sep="\t")

#These stats can be used on the command prompt to check data
#FCstats <- boxplot.stats(Averaged_Orderd_KO_data[,1])
#FCstats$conf

####################################################################
####################################################################
####################################################################

#SECTION 3
#subset data to sample only PRC2-target genes

#read in geneIDs of PRC2 target genes
Prc2targets <- read.table("./JGI_ListofK27SilentGenesNCUs.txt", header=FALSE,stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

###SUBSET DATA FOR ALL KO SAMPLES
Prc2targetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%Prc2targets[,1])

###SUBSET DATA FOR AVERAGED KO SAMPLES
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%Prc2targets[,1])

###calculateSummaryStats for PCR2 target gene data
TPMmean <- colMeans(AVERAGE_Prc2targetTPM)
TPMmedian <- colMedians(AVERAGE_Prc2targetTPM)
TPMstdev <- apply(AVERAGE_Prc2targetTPM,2,sd)
w<- cbind(TPMmean,TPMmedian,TPMstdev)

#write.table(w, file="~/Dropbox/DropBOX Documents/Zack Papers/2018_RemodellerPaper/IswiPaperFilesFromMasayuki/SummaryTablesFromGenomics/KO_AverageTPM_stats.txt", sep="\t")

#these commands can be used pn the command line to check stats of data
#FCstats <- boxplot.stats(AVERAGE_Prc2targetTPM[,12])
#FCstats$conf

#melt data to get it into a format ggplots can use
library(reshape2)

meltedPRC2targetData <- melt(Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

#################################################################################

##PLOT AVERAGE DATA FOR FIGURE 1A

##################################################################
#melt the data and label value "Count" and the sample IDs GeneID and Sample
#meltedData <- melt(tpm, value.name = "Count",
#varnames=c('GeneID', 'Sample'))

library(reshape2)
meltedAveragePRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

altorder = rev(c( "WT","set-7","cac-1","cac-2","cac-3","naf-1","naf-2","asf-1","ATRX"))
meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample, altorder)

# this works for all genes but formatting is terrible
library(ggplot2)
xlabels = averageRowIDs
colors = rev(c( "#4575b4","#fee090","#fee090", "#fee090", "#fee090", "#4575b4", "#4575b4", "#4575b4", "#4575b4"))

box<-ggplot(meltedAveragePRC2targetData, aes(x=Sample, y=Count)) +
  labs(y="Expression Level (Transcripts per Million)", x="Strain") +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill=colors, size=0.1, coef=1.5, lwd=5) +
  coord_flip(ylim=c(-1,25))+
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

ggsave(filename = "histone_chaperone_boxplot_rerun.pdf", plot = box, dpi=600, height= 3, width=4)


########################3
############################################################################

library(pheatmap)
library(RColorBrewer)

breaks1=seq(-4, 5, by=.09) #This is to set a custom heatmaps scale. Not used here.

#filter all data with no changes >>>>> if the sum of all rows divided by the # of rows is > 0.

##Experiment with different heatmap strategies
##strategy 1 - Plot the average TMP of all genes that change expression; row normalization; clustredRows

##Keep Only Genes that are expressed in at least one sample
GenesWithChanges <- subset(allgene_fac, (rowSums(fac_TPM) > 0))
GenesWithChanges_95 <- subset(allgene_chromatin_reg, (rowSums(Prc2targetTPM) > 0))

#plot in the desired column order; did this by subsetting the dataset based on sample list 'altorder' above

heatmap<- pheatmap(GenesWithChanges[,rev(altorder)], color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                   cellwidth = NA, cellheight = NA,  cluster_rows =T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                   legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot<-heatmap[[4]]

ggsave(filename = "./fac_het_tpm.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)
#dev.off()
#clustering_method="centroid", clustering_distance_cols="euclidean", 

#######################################################################
#######################################################################

allgene_chromatin_reg = subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data) %in% chromatin_regulators$Gene.ID)

meltedAverage_chromatin_reg <- melt(allgene_chromatin_reg, value.name = 'Count',
                                    varnames=c('GeneID', 'Sample'))

meltedAverage_chromatin_reg$Sample <- factor(meltedAverage_chromatin_reg$Sample, altorder)

chromatin_reg_TPM <- subset(allDataTPM, rownames(allDataTPM)%in%chromatin_regulators$Gene.ID)

GenesWithChanges = GenesWithChanges + 1
log_GenesWithChanges = log2(GenesWithChanges)

#####
allgene_fac = subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data) %in% fac_het$Gene.ID)

meltedAverage_fac <- melt(allgene_fac, value.name = 'Count',
                                    varnames=c('GeneID', 'Sample'))

meltedAverage_fac$Sample <- factor(meltedAverage_fac$Sample, altorder)

fac_TPM <- subset(allDataTPM, rownames(allDataTPM)%in%fac_het$Gene.ID)

GenesWithChanges = GenesWithChanges + 1
log_GenesWithChanges = log2(GenesWithChanges)


draw.triple.venn(area1 = 1711, area2 = 1625, area3 = 1953, n12 = 1261, n23 = 947, n13 = 1016, n123 = 823,
                 category = c("cac-1", "cac-2", "cac-3"), fill = c("Blue", "Red", "Green"))

