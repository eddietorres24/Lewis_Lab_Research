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
library("scater")
library(tidyr)

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
countdataInteractors <- read.table("./text_files/readcounts_All_CAF1paper.merged.txt", header=TRUE,stringsAsFactors=FALSE, row.names=1, check.names=FALSE, sep="\t")
###MAY NEED TO SKIP ONE ROW w/ skip =1 

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

#use to find correct columns
srrWT  <- c("SRR8269825","SRR8269775","SRR8269782","SRR8269810")
srrcac1  <- c("149-115_RNA_cac-1__Rep1_S139","149-116_RNA_cac-1__Rep2_S140","149-117_RNA_cac-1__Rep3_S141")
srrcac2  <- c("SRR7970629","SRR7970630","SRR7970631")
srrset7  <- c("SRR10916163","SRR10916164","SRR10916165")
srrash1  <- c("SRR7690267","SRR7690268")

cn  <- colnames(allDataTPM)
bn  <- basename(cn)
ids <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", bn)

idx <- which(ids %in% srrash1)   # <- column numbers you want
idx
bn[idx]                        # readable column names (optional)
# subset if needed:
sub <- allDataTPM[, idx, drop = FALSE]

#move the subset of genes you want to plot into a new matrix
sampleNames <- colnames(allDataTPM)
write.table(sampleNames, file="names.txt")
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,c(6,8,13,19)], na.rm = TRUE),
                                 rowMeans(allDataTPM[,c(34,36,39)], na.rm = TRUE),
                                 rowMeans(allDataTPM[,c(7,10,23)], na.rm = TRUE),
                                 rowMeans(allDataTPM[,c(1,2)], na.rm = TRUE))
averageRowIDs=c("WT","cac-1","set-7","ash-1")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

####################################################################

#SECTION 3
#subset data to sample only PRC2-target genes (or other gene sets of interest)

#read in geneIDs of PRC2 target genes
#reading in csvs w/ upregulated genes in CAF-1 mutants
Prc2targets <- read.table("./bed_files/H3K27me3_methylated_genes_FINAL.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

###SUBSET DATA FOR ALL KO SAMPLES
Prc2targetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%Prc2targets[,10])
#genes i had to rename in my bed file due to spaces
hsp <- allDataTPM[rownames(allDataTPM) == "hsp30-like protein", ]
lac <- allDataTPM[rownames(allDataTPM) == "laccase precursor", ]
#add back missing genes
Prc2targetTPM <- rbind(Prc2targetTPM, hsp, lac)
#rename rows of added genes
rownames(Prc2targetTPM)[rownames(Prc2targetTPM) == "hsp"] <- "hsp30-like protein"
rownames(Prc2targetTPM)[rownames(Prc2targetTPM) == "lac"] <- "laccase precursor"

###SUBSET DATA FOR AVERAGED KO SAMPLES
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%rownames(Prc2targetTPM))

###Subset data to filter out non-repressed PRC2 targets regions that got through (cutting out any genes over 10 tpm in WT).
###These are likely genes whose promoter are not marked by H3K27me3, genes on the edge of K27 regions, and/or bivalent genes
AVERAGE_Prc2targetTPM <- subset(AVERAGE_Prc2targetTPM, (AVERAGE_Prc2targetTPM[,1] < 2))
AVERAGE_AlldataTPM <- subset(Averaged_Orderd_KO_data, (Averaged_Orderd_KO_data[,1] > -0.1))

###resubet PRC2 targets after filtering
Prc2targetTPM <- subset(Prc2targetTPM, rownames(Prc2targetTPM)%in%rownames(AVERAGE_Prc2targetTPM))

###Add sudocount and log transform (if necessary)
AVERAGE_Prc2targetTPM <- AVERAGE_Prc2targetTPM + 1
AVERAGE_Prc2targetTPM <- log2(AVERAGE_Prc2targetTPM)

AVERAGE_AlldataTPM <- AVERAGE_AlldataTPM + 1
AVERAGE_AlldataTPM <- log2(AVERAGE_AlldataTPM)
#melt data to get it into a format ggplots can use (load library "reshape2")
library(reshape2)

meltedAveragePRC2targetData <- melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                             varnames=c('GeneID', 'Sample'))

#Alldata
meltedAverageAllData <- melt(AVERAGE_AlldataTPM, value.name = 'Count',
                      varnames=c('GeneID', 'Sample'))



#################################################################################

##PLOT AVERAGE DATA FOR FIGURE 1A

##################################################################
#melt the data and label value "Count" and the sample IDs GeneID and Sample
#meltedData <- melt(tpm, value.name = "Count",
#varnames=c('GeneID', 'Sample'))

altorder = rev(c("WT","cac1","set7","ash1"))

meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample)

meltedAverageAllData$Sample <- factor(meltedAverageAllData$Sample)

# Plot box & whisker chart
library(ggplot2)
xlabels = averageRowIDs
colors = c( "#4575b4","#fee090","#fee090", "#fee090", "#fee090")

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

############################################

###Make histogram to visualize distribution of data

df <- as.data.frame(AVERAGE_Prc2targetTPM)

long_df <- df %>%
  select(WT, `cac-1`) %>%
  pivot_longer(cols = everything(), names_to = "Strain", values_to = "TPM")

# Plot overlapping histograms
ggplot(long_df, aes(x = TPM, fill = Strain)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, color = "black") +
  theme_minimal() +
  labs(title = "TPM Distribution in WT and cac-1",
       x = "TPM",
       y = "Frequency") +
  scale_fill_manual(values = c("WT" = "steelblue", "cac-1" = "tomato"))

#Histogram shows non-normal distribution, so proceed with wilcox test
#t-test for significance in difference of mean expression values of PRC2 targets b/w WT & mutants

set7_stat <- wilcox.test(df$`set-7`, df$WT, alternative = "greater")
cac1_stat <- wilcox.test(df$`cac-1`, df$WT, alternative = "greater")
cac2_stat <- wilcox.test(df$`cac-2`, df$WT, alternative = "greater")
# cac1_2_tstat <- wilcox.test(df$`cac-1_2`, df$WT, alternative = "greater")
cac3_stat <- wilcox.test(df$`cac-3`, df$WT, alternative = "greater")
naf1_stat <- wilcox.test(df$`naf-1`, df$WT, alternative = "greater")
naf2_stat <- wilcox.test(df$`naf-2`, df$WT, alternative = "greater")
asf1_stat <- wilcox.test(df$`asf-1`, df$WT, alternative = "greater")
atrx_stat <- wilcox.test(df$`ATRX`, df$WT, alternative = "greater")

###Make t-test dataframe

#Make a list of t-test results objects
wilcox_list <- list(set7_stat, cac1_stat, cac2_stat, cac3_stat, naf1_stat, naf2_stat, asf1_stat, atrx_stat)
names(wilcox_list) <- c("set7", "cac1", "cac2", "cac3", "naf1", "naf2", "asf1", "atrx")

# Function to extract relevant info from each t-test
# Name each element in the list according to the strain
names(wilcox_list) <- c("set-7", "cac-1", "cac-2", "cac-3", "naf-1", "naf-2", "asf-1", "ATRX")

# Extract relevant stats into a dataframe
wilcox_df <- do.call(rbind, lapply(wilcox_list, function(res) {
  data.frame(
    statistic = unname(res$statistic),
    p_value = res$p.value,
    alternative = res$alternative,
    method = res$method
  )
}))

# Set strain names as row names
rownames(wilcox_df) <- names(wilcox_list)

# Save to CSV
write.csv(wilcox_df, "Wilcoxon_summary_stats.csv", row.names = TRUE)

################################################

#I want to make a violin plot instead of boxplot. Going to adapt Abby's code
library(grDevices)

# Define factor order and formatted x-axis labels
label_order <- c("WT", "cac-1", "cac-2", "set-7", "ash-1")
formatted_labels <- c(
  "WT",
  expression(italic("\u0394cac-1")),
  expression(italic("\u0394cac-2")),
  expression(italic("\u0394set-7")),
  expression(italic("\u0394ash-1"))
)

# Group and summarize
total2 <- meltedAveragePRC2targetData %>%
  group_by(Sample)

total_dist <- meltedAveragePRC2targetData %>%
  group_by(Sample) %>%
  summarise(num = n())

# Set dodge for violin spacing
dodge <- position_dodge(width = 1)

# Build plot
violin <- total2 %>%
  left_join(total_dist) %>%
  arrange(factor(Sample, levels = label_order)) %>%
  mutate(Sample = factor(Sample, levels = label_order)) %>%
  ggplot(aes(x = Sample, y = Count)) +
  geom_violin(position = dodge, scale = "width", trim = FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
  scale_x_discrete(labels = formatted_labels) +
  labs(
    x = "Strain",
    y = expression("Expression Level (log"[2]~"(TPM+1))")) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black")
  )

print(violin)

ggsave("./Histone_Chap_PRC2_Gene_expression_FINAL_TEST.pdf", plot=violin, width = 10, height = 8, unit="in",  dpi=400)

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

heatmap <- pheatmap(K27_not_lost, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                    cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                    legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#to plot with ggplot, you need to extract [[4]] from the heatmap object
heatmap_plot <- heatmap[[4]]

###############################################
#replicate scaling for other heatmaps
library(pheatmap)   # for the internal helper

rna_cols <- c("WT","cac-1","set-7","ash-1")
M <- as.matrix(GenesWithChanges[, rna_cols, drop = FALSE])

# identical to pheatmap(scale="row")
rna_z <- pheatmap:::scale_rows(M)

# (optional) clip for display like pheatmap often does
rna_z <- pmin(pmax(rna_z, -3), 3)
###############################################

ggsave(filename = "./CAF-1_K27_Paper_Final.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)
#dev.off()
