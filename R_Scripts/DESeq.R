---
title: "DESeq2"
author: "Abby Deaven"
date: "2022-10-30"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "C:/Users/abbyd/OneDrive - University of Georgia/Lewis Lab/Research/Data/Development/Bioinformatics/RNAseq/RNAseqPlotting_New10-2022/res")
```a


## R Markdown
This code generates log2foldchange expression data from RNAseq data that have been mapped using the featurecounts command in the Subread package.

```{r packages}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
BiocManager::install("DESeq2")
BiocManager::install("scales")
BiocManager::install("IHW")

library("DESeq2")
library("ggplot2")
library("scales")
library(readr)
library(tidyr)
library(tidyverse)
library(dplyr)


## Importing featurecounts matrix
Read in the featurecounts output, convert to a matrix, and remove unneccesary columns to include ONLY the count numbers, then convert columns to numeric.

setwd("C:\\Users\\eddie\\Desktop\\Research\\GitHub")
cts <- as.matrix(read.table("./text_files/readcounts_FINAL.txt", sep="\t", row.names =1, header = TRUE))
cts <- as.matrix(cts[,-c(1:5)])
cts_r <- as.matrix(read.table("./text_files/readcounts_rtt109.txt", sep="\t", row.names =1, header = TRUE))
cts_r <- as.matrix(cts_r[,-c(1:5)])
cts_f <- cbind(cts, cts_r)
## change ncol to match number of columns
cts_numeric <- matrix(as.numeric(cts_f), ncol = 31)
dimnames(cts_numeric) <- list(rownames(cts_f), colnames(cts_f))




#Export column names to a text file, open in a text editor and replace SRRs with short/identifiable names, save, and read in the file again to replace column titles with human readable names.


samplesname <- colnames(cts_numeric)
write.table(samplesname, file="samplenames.txt", sep="\t")

#after editing sample names in excel or text editor (new names are in column 3), read in the spreadsheet with hand annotated sample names that are short 
rename <- read.table("./samplenames.txt", header=FALSE, skip=1, sep = "\t")
#assign short sample names to variable
cnames <- rename$V2
#replace column names
colnames(cts_numeric) <- cnames


#All samples must have the same number of replicates. I'm averaging across samples that have multiple technical replicates, then reducing the number of replicates to 2 for each condition.


cac1 <- rowMeans(cts_numeric[, 1:3], na.rm = TRUE)
cac2 <- rowMeans(cts_numeric[, 4:6], na.rm = TRUE)
cac3 <- rowMeans(cts_numeric[, 7:9], na.rm = TRUE)
WT <- rowMeans(cts_numeric[, 22:25], na.rm = TRUE)
set7 <- rowMeans(cts_numeric[, 26:28], na.rm = TRUE)
rtt109 <- rowMeans(cts_numeric[, 29:31], na.rm = TRUE)

## Importing DEseq metadata
#The coldata file contains the column names that are shown in the cts matrix, the sample identifier (for pooling replicates), and whether each sample is a control or experimental variable. For example:

# head(coldata)
#                X   condition         type
#1        Hy2489_2 Hyphae_48hr experimental
#2        Hy2489_1 Hyphae_48hr experimental
#3        wt_rep_1   WTMycelia      control
#4        wt_rep_2   WTMycelia      control
#5 NoCarbon_1_Rep1    NoCarbon experimental
#6 NoCarbon_1_Rep2    NoCarbon experimental

#After that, convert each column to a factor.


coldata2 <- read.csv("./csv_files/coldata_caf.csv", header= TRUE, row.names = 1)
coldata2$condition <- factor(coldata2$condition)
coldata2$type <- factor(coldata2$type)


#Confirm that the number and order of samples are the same in each dataset before proceeding.

all(rownames(coldata2) %in% colnames(cts_numeric))
all(rownames(coldata2) == colnames(cts_numeric))

cts_numeric_3 = cts_numeric + 1

##### CREATE DESEQ DATASET AND RUN DESEQ2 #####

#Before running DESeq2, create a matrix from the objects you created and specify the levels you'll be comparing (in this case, experimental vs. control type, separated by each individual condition. Then, pre-filter to remove very low read genes.)


##create DEseq dataset
dds2<- DESeqDataSetFromMatrix(countData = round(cts_numeric_3),
                              colData = coldata2,
                              design = ~ condition)
dds2

##pre-filteringreads <10
keep <- rowSums(counts(dds2)) > 19
dds2<- dds2[keep,]

#specify level of comparison
dds2$condition <- relevel(dds2$condition, ref = "WT")
```

#Now, run DEseq and plot dispesion estimates to see the quality of the fit. 


dds2 <- DESeq(dds2)
plotDispEsts(dds2)

#Export genes with statistically significant differential expression when compared against the wild type.


alpha = 0.05

cac1 <- results(dds2, alpha=alpha, contrast=c("condition", "cac1", "WT"))
cac2 <- results(dds2, alpha=alpha, contrast=c("condition", "cac2", "WT"))
cac3 <- results(dds2, alpha=alpha, contrast=c("condition", "cac3", "WT"))
#naf1 <- results(dds2, alpha=alpha, contrast=c("condition", "naf1", "WT"))
#naf2 <- results(dds2, alpha=alpha, contrast=c("condition", "naf2", "WT"))
#asf1 <- results(dds2, alpha=alpha, contrast=c("condition", "asf1", "WT"))
#ATRX <- results(dds2, alpha=alpha, contrast=c("condition", "ATRX", "WT"))
set7 <- results(dds2, alpha=alpha, contrast=c("condition", "set7", "WT"))
rtt109 <- results(dds2, alpha=alpha, contrast=c("condition", "rtt109", "WT"))


#Subset to only export significant differentially expressed genes

cac1 <- subset(cac1, padj < 0.05)
cac2 <- subset(cac2, padj < 0.05)
cac3 <- subset(cac3, padj < 0.05)
#naf1 <- subset(naf1, padj < 0.05)
#naf2 <- subset(naf2, padj < 0.05)
#asf1 <- subset(asf1, padj < 0.05)
#ATRX <- subset(ATRX, padj < 0.05)
set7 <- subset(set7, padj < 0.05)
rtt109 <- subset(rtt109, padj < 0.05)

write.csv(as.data.frame(cac1), file = "cac1.csv")
write.csv(as.data.frame(cac2), file = "cac2.csv")
write.csv(as.data.frame(cac3), file = "cac3.csv")
#write.csv(as.data.frame(naf1), file = "naf1.csv")
#write.csv(as.data.frame(naf2), file = "naf2.csv")
#write.csv(as.data.frame(asf1), file = "asf1.csv")
#write.csv(as.data.frame(ATRX), file = "ATRX.csv")
write.csv(as.data.frame(set7), file = "set7.csv")
write.csv(as.data.frame(rtt109), file = "rtt109.csv")

### Analyzing DE Genes

#Read results files back into r as a single matrix


knitr::opts_knit$set(root.dir = "C:\\Users\\eddie\\Desktop\\Research\\GitHub\\Research")

setwd("C:\\Users\\eddie\\Desktop\\Research\\GitHub\\Research")
list_of_files <- list.files(path = "./cac_DEseq",
                            recursive = TRUE,
                            pattern = ".csv$")

induced_genes <- readr::read_csv(list_of_files, id = "file_name")

##convert to wide format
l2fc_data <- data.frame(pivot_wider(data = induced_genes, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))
rownames(l2fc_data) <- l2fc_data[,1]
l2fc_data <- l2fc_data[,-1]
l2fc_data[is.na(l2fc_data)] <- 0
#rownames(l2fc_data) <- l2fc_data$gene


library(tibble)
l2fc <- tibble::rownames_to_column(l2fc_data, "VALUE")


#Find genes upregulated at any stage of perithecial development


peri <- l2fc_data[,2:5]
peri_induced <- peri  %>% filter_at(vars(1:4), any_vars(. >4))
peri_3 <- peri %>% filter_at(vars(1), any_vars(. >4))
peri_4 <- peri %>% filter_at(vars(2), any_vars(. >4))
peri_5 <- peri %>% filter_at(vars(3), any_vars(. >4))
peri_6 <- peri %>% filter_at(vars(4), any_vars(. >4))

write.csv(as.data.frame(peri_induced), file = "peri_induced.csv")


#Compare gene upregulation using UpSet plots. Data is entered by first listing the categories and the number of individuals in each category. After that, list intersections. I like to do this by calculating the number of overlaps directly in the input, but you can also do this separately. For example, if X has 10 samples, Y has 15 samples, and Z has 20 samples, you would enter like this:

#input <- c(
#  X = 10,
#  Y = 15,
#  Z = 20,
#  "X&Y"= nrow(intersect(X, Y)),  "X&Z"= nrow(intersect(X, Z)),
#  "Y&X"= nrow(intersect(Y, Z)),
#  "X&Y&Z" = nrow(intersect(intersect(X, Y),Z)))
#)


#To make the UpSet plot, run the code "upset(fromExpression(input)," and input parameters (read the documentation for full list of options). You MUST specify the number of categories (nsets) and number of intersections (nintersects), which INCLUDES the number of individual categories.


install.packages("UpSetR")
library(UpSetR)

library(UpSetR)

# Dataset
input <- c(
  Day3 = 1233,
  Day4 = 1407,
  Day5 = 1457,
  Day6 = 1489,
  "Day3&Day4" = nrow(intersect(peri_3, peri_4)),
  "Day3&Day5" = nrow(intersect(peri_3, peri_5)),
  "Day3&Day6" = nrow(intersect(peri_3, peri_6)),
  "Day4&Day5" = nrow(intersect(peri_4, peri_5)),
  "Day4&Day6" = nrow(intersect(peri_4, peri_6)),
  "Day5&Day6" = nrow(intersect(peri_5, peri_6)),
  "Day3&Day4&Day5" = nrow(intersect(intersect(peri_3, peri_4), peri_5)),
  "Day4&Day5&Day6" = nrow(intersect(intersect(peri_4, peri_5), peri_6)),
  "Day3&Day5&Day6" = nrow(intersect(intersect(peri_3, peri_5), peri_6)),
  "Day3&Day4&Day6" = nrow(intersect(intersect(peri_3, peri_4), peri_6)),
  "Day3&Day4&Day5&Day6" = nrow(intersect(intersect(intersect(peri_3, peri_4), peri_6), peri_5)))
  
 upset(fromExpression(input), 
      nintersects = 15, 
      nsets = 4, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.4, 
      point.size = 2.8, 
      line.size = 1
      )

#IGV Track

ngenes = read.csv("neurospora_genes.csv")

combine <- transform(merge(ngenes, l2fc_data, by = "VALUE"))
track = combine %>% select("SequenceID", "FeatureStart", "FeatureEnd", "VALUE", "asf1.csv", "ATRX.csv", "cac1.csv", "cac2.csv", "cac3.csv", "naf1.csv", "naf2.csv", "set7.csv")

write.table(track, file="cactrack.igv", sep="\t", row.names = FALSE, quote = FALSE)


