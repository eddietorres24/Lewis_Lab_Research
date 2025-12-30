#### Violin plot script ###


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("grDevices")

library(dplyr)
library(ggplot2)
library(grDevices)


#### Pivoting matrix of 3dpf, 4dpf, 5dpf, 6dpf, and mycelia samples with log2(TPM+1) longer ###
PRC2violin <- AVERAGE_Prc2targetTPM %>% pivot_longer(names_to = "Sample", values_to = "logTPM_plusPC")

### Subsetting to H3K27me3 and non-H3K27me3 genes (dev_non_prc2 = list of non-H3K27me3 marked but developmentally-induced genes) and then adding a column stating they are H3K27me3-depleted
dev_non_prc2_longer <- data.frame(dev_samples_log_longer[dev_samples_log_longer$GeneID %in% dev_non_prc2$x,])
dev_non_prc2_longer$Category <- "H3K27me3-Depleted"

#Same as above but for H3K27me3-enriched developmentally induced genes
prc2_longer <- dev_samples_log_longer[dev_samples_log_longer$GeneID %in% prc2_induced$x,]
prc2_longer$Category <- "H3K27me3-Enriched"

#merging each matrix together
total <- rbind(dev_non_prc2_longer, prc2_longer)

### this part is just for determining the order of samples in the violin plot -- groups genes by sample (3dpf, mycelia, etc) and then summarizes how many total genes are in each sample (is the same for each one) ###
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
  arrange(factor(Sample, levels = c("WT", "set-7", "cac-1", "cac-2", "cac-3", "naf-1", "naf-2", "asf-1", "ATRX"))) %>%
  mutate(Sample = factor(Sample)) %>%
  ### everything below is the actual violin plot ###
  ggplot(aes(x=Sample, y=Count)) + 
  geom_violin(position = dodge, scale="width", trim=FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  scale_fill_manual("hello",values = c("orchid1", "springgreen3")) +
  labs(x = "Strain",y = expression("Expression Level (log "[2]~"(TPM+1)")) + 
  theme_classic(base_size = 20)

print(violin)
ggsave("./Histone_Chape_PRC2_Gene_expression.pdf", plot=violin, width = 10, height = 8, unit="in",  dpi=400)
