library(ggplot2)
library(tidyverse)
library(fs)
library(readxl)
library(dplyr)
library(expss)
library(RVAideMemoire)
library(ggprism)
library(ggrepel)
library(WGCNA)
library(GO.db)
library(genefilter)
library(sparseMatrixStats)
library(DelayedMatrixStats)
library(tximport)
library(tximportData)
library(DESeq2)
library(rhdf5)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(data.table)
library(ashr)

#setwd("~/Downloads/data")
####Cryptic species and mortality analysis####
#Scatterplot of cryptic species designations

cryptic <- read_excel("./crypticdata.xlsx", na = "")
ggplot(cryptic, aes(x=PC1, y=PC2, color=mortality, shape=CrypticSpeciesLabel, label=ID)) +
  geom_vline(aes(xintercept=0.03),linetype="dashed",color="black") +
  geom_vline(aes(xintercept=-0.020),linetype="dashed",color="black") +
  geom_vline(aes(xintercept=-0.080),linetype="dashed",color="black") +
  geom_vline(aes(xintercept=-0.1175),linetype="dashed",color="black") +
  geom_point(size=5) +
  scale_shape_manual(values=c(15,17,6,18,19)) +
  scale_color_manual(values=c("seagreen4","orange4","gray30")) +
  geom_text_repel(nudge_x=0.005, nudge_y=0.005,max.overlaps=80,size=6) + 
  labs(x = "PC 1",
       y = expression("PC 2")) +
  theme_bw() + theme(legend.position = "none", axis.text = element_text(size=20), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_text(size=20), panel.border = element_blank(), axis.line = element_line(colour = "black"))

#map of the study site 
#with cryptic species and relative heat stress resistance designation
library(sf)

data <- read_csv("./latlongfile.csv")
points = st_as_sf(data, coords = c("X","Y"), crs = 4326)
plot(st_geometry(points), pch=16,col="navy")

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale="medium", returnclass="sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(134.5,134.515), ylim = c(7.285,7.3), expand = FALSE) +
  theme_bw()

##Combined points and coordinates map with reef cutouts in Illustrator.

##Fisher's exact test for geography
csp_map <- matrix(c(17,14,2,5),ncol=2,dimnames=list(c("Csp5", "Csp2"),c("PR7","PR9")),byrow=TRUE)
fisher.test(csp_map)

##Fisher's exact test for mortality
csp_analysis <- matrix(c(26,5,1,6),ncol=2,dimnames=list(c("Csp5", "Csp2"),c("Alive","Dead")),byrow=TRUE)
fisher.test(csp_analysis)

csp_plot <- matrix(c(26,5,1,6),ncol=2)
colnames(csp_plot) <- c("Csp5","Csp2")
rownames(csp_plot) <- c("Alive", "Dead")
library("reshape2")
tmp <- melt(csp_plot)
names(tmp) <- c("Mortality", "CrypticSpecies", "NumberGenets")

ggplot(tmp, aes(x=CrypticSpecies, y=NumberGenets, fill=Mortality)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  scale_fill_manual(values=c("seagreen4","orange4")) +
  theme_bw()

####Data set 1:pre heat stress ALL####
dir <- ("~/Downloads/data/TranscriptomeData")
list.files(dir)
samples_pre <- read_excel(file.path(dir, "sampleIDs_data.xlsx"), 1)
files_pre <- file.path(dir, "Samples_Current", samples_pre$SampleID, "abundance.h5")
names(files_pre) <- paste0("sample", 1:58)
all(file.exists(files_pre))
txi_pre <- tximport(files_pre, type = "kallisto", txOut = TRUE)
head(txi_pre$counts)

sampleTable_pre <- data.frame(samples_pre$SampleID, samples_pre$HR_Cat, samples_pre$Cryptic, samples_pre$AverageDaysinHS)
rownames(sampleTable_pre) <- colnames(txi_pre$counts)
all(rownames(sampleTable_pre) %in% colnames(txi_pre$counts))
all(rownames(sampleTable_pre) == colnames(txi_pre$counts))
colnames(sampleTable_pre) <- c("SampleID", "HR_Cat", "Cryptic", "AverageDaysinHS")

sampleTable_pre$Cryptic <- as.factor(sampleTable_pre$Cryptic)
dds_pre <- DESeqDataSetFromTximport(txi_pre, colData = sampleTable_pre, design = ~ HR_Cat + Cryptic)
dds_pre$HR_Cat <- factor(dds_pre$HR_Cat, levels=c("Dead","Low","Moderate","High"))
dds_pre <- DESeq(dds_pre)

resultsNames(dds_pre)
####Different pairwise results pre by cryptic species####
res_cryptic_pre <- results(dds_pre, contrast=c("Cryptic", "2", "5"),
                           independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_cryptic_pre
res_cryptic_pre <- res_cryptic_pre[order(res_cryptic_pre$pvalue),]
#Create a Folder Titled R_Files
write.csv(res_cryptic_pre,"./R_Files/ALL_cryptic_pre_raw.csv")
summary(res_cryptic_pre)
metadata(res_cryptic_pre)$filterThreshold
sum(res_cryptic_pre$padj < 0.05, na.rm = TRUE)

res_Crp_pre_LFC <- lfcShrink(dds_pre, contrast=c("Cryptic", "2", "5"), res=res_cryptic_pre, type="ashr")
df_Crp_pre_LFC <- as.data.frame(res_Crp_pre_LFC)
df_Crp_pre_LFC <- df_Crp_pre_LFC[order(df_Crp_pre_LFC$pvalue), ]
summary(res_Crp_pre_LFC, alpha = 0.05)
sum(df_Crp_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_Crp_pre_LFC, "./R_Files/ALL_cryptic_pre_ADJUSTED.csv")

####Different pairwise results pre####
res_LM_pre <- results(dds_pre, contrast=c("HR_Cat", "Low", "Moderate"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_LM_pre
res_LM_pre <- res_LM_pre[order(res_LM_pre$pvalue),]
write.csv(res_LM_pre,"./R_Files/ALL_LM_pre_raw.csv")
summary(res_LM_pre)
sum(res_LM_pre$padj < 0.05, na.rm = TRUE)
#73 out of 27964
#LFC > 0 29
#LFC < 0 44
res_LM_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "Low", "Moderate"), res=res_LM_pre, type="ashr")
df_LM_pre_LFC <- as.data.frame(res_LM_pre_LFC)
df_LM_pre_LFC <- df_LM_pre_LFC[order(df_LM_pre_LFC$pvalue), ]
summary(res_LM_pre_LFC, alpha = 0.05)
sum(df_LM_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_LM_pre_LFC, "./R_Files/ALL_LM_pre_ADJUSTED.csv")

res_LH_pre <- results(dds_pre, contrast=c("HR_Cat", "Low", "High"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_LH_pre
res_LH_pre <- res_LH_pre[order(res_LH_pre$pvalue),]
write.csv(res_LH_pre,"./R_Files/ALL_LH_pre_raw.csv")
summary(res_LH_pre)
sum(res_LH_pre$padj < 0.05, na.rm = TRUE)
#1558 out of 27964
#LFC > 0 420
#LFC < 0 1138
res_LH_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "Low", "High"), res=res_LH_pre, type="ashr")
df_LH_pre_LFC <- as.data.frame(res_LH_pre_LFC)
df_LH_pre_LFC <- df_LH_pre_LFC[order(df_LH_pre_LFC$pvalue), ]
summary(res_LH_pre_LFC, alpha = 0.05)
sum(df_LH_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_LH_pre_LFC, "./R_Files/ALL_LH_pre_ADJUSTED.csv")


res_LD_pre <- results(dds_pre, contrast=c("HR_Cat", "Low", "Dead"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_LD_pre
res_LD_pre <- res_LD_pre[order(res_LD_pre$pvalue),]
write.csv(res_LD_pre,"./R_Files/ALL_LD_pre_raw.csv")
summary(res_LD_pre)
sum(res_LD_pre$padj < 0.05, na.rm = TRUE)
#57 out of 27964
#LFC > 0 37
#LFC < 0 20
res_LD_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "Low", "Dead"), res=res_LD_pre, type="ashr")
df_LD_pre_LFC <- as.data.frame(res_LD_pre_LFC)
df_LD_pre_LFC <- df_LD_pre_LFC[order(df_LD_pre_LFC$pvalue), ]
summary(res_LD_pre_LFC, alpha = 0.05)
sum(df_LD_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_LD_pre_LFC, "./R_Files/ALL_LD_pre_ADJUSTED.csv")

res_MH_pre <- results(dds_pre, contrast=c("HR_Cat", "Moderate", "High"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_MH_pre
res_MH_pre <- res_MH_pre[order(res_MH_pre$pvalue),]
write.csv(res_MH_pre,"./R_Files/ALL_MH_pre_raw.csv")
summary(res_MH_pre)
sum(res_MH_pre$padj < 0.05, na.rm = TRUE)
#24 out of 27964
#LFC > 0 7
#LFC < 0 17
res_MH_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "Moderate", "High"), res=res_MH_pre, type="ashr")
df_MH_pre_LFC <- as.data.frame(res_MH_pre_LFC)
df_MH_pre_LFC <- df_MH_pre_LFC[order(df_MH_pre_LFC$pvalue), ]
summary(res_MH_pre_LFC, alpha = 0.05)
sum(df_MH_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_MH_pre_LFC, "./R_Files/ALL_MH_pre_ADJUSTED.csv")

res_MD_pre <- results(dds_pre, contrast=c("HR_Cat", "Moderate", "Dead"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_MD_pre
res_MD_pre <- res_MD_pre[order(res_MD_pre$pvalue),]
write.csv(res_MD_pre,"./R_Files/ALL_MD_pre_raw.csv")
summary(res_MD_pre)
sum(res_MD_pre$padj < 0.05, na.rm = TRUE)
#58 out of 27964
#LFC > 0 51
#LFC < 0 7
res_MD_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "Moderate", "Dead"), res=res_MD_pre, type="ashr")
df_MD_pre_LFC <- as.data.frame(res_MD_pre_LFC)
df_MD_pre_LFC <- df_MD_pre_LFC[order(df_MD_pre_LFC$pvalue), ]
summary(res_MD_pre_LFC, alpha = 0.05)
sum(df_MD_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_MD_pre_LFC, "./R_Files/ALL_MD_pre_ADJUSTED.csv")

res_HD_pre <- results(dds_pre, contrast=c("HR_Cat", "High", "Dead"),
                      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
res_HD_pre
res_HD_pre <- res_HD_pre[order(res_HD_pre$pvalue),]
write.csv(res_HD_pre,"./R_Files/ALL_HD_pre_raw.csv")
summary(res_HD_pre)
sum(res_HD_pre$padj < 0.05, na.rm = TRUE)
#110 out of 27964
#LFC > 0 106
#LFC < 0 4
res_HD_pre_LFC <- lfcShrink(dds_pre, contrast=c("HR_Cat", "High", "Dead"), res=res_HD_pre, type="ashr")
df_HD_pre_LFC <- as.data.frame(res_HD_pre_LFC)
df_HD_pre_LFC <- df_HD_pre_LFC[order(df_HD_pre_LFC$pvalue), ]
summary(res_HD_pre_LFC, alpha = 0.05)
sum(df_HD_pre_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_HD_pre_LFC, "./R_Files/ALL_HD_pre_ADJUSTED.csv")

####Fig S1: Supp PRE Volcano Plots####
#Make a Figures folder
pdf("./Figures/S1_Baseline_Expression_Comparisons.pdf", width = 11, height = 8.5)
par(mfrow=c(2,3), mar=c(4,4,2,1), oma = c(0, 0, 3, 0))
xlim <- c(1,1e4); ylim <- c(-3,3); 
plotMA(res_LM_pre_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Low vs. Moderate")
plotMA(res_LD_pre_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Low vs. Dead")
plotMA(res_MH_pre_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Moderate vs. High")
plotMA(res_MD_pre_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Moderate vs. Dead")
plotMA(res_HD_pre_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "High vs. Dead")

mtext("Baseline Expression Comparisons Between Heat Resistance Categories", 
      outer = TRUE, 
      side = 3, 
      cex = 1.5, 
      font = 2, 
      line = 1)

dev.off()

####Main Fig 2: PCA + Volcano Plots PRE####
window_xlimit1 <- 4.5
window_ylimit1 <- 4
df_M2A <- as.data.frame(res_Crp_pre_LFC)
df_M2A$Direction <- "Not Significant"
df_M2A$Direction[df_M2A$padj < 0.05 & df_M2A$log2FoldChange > 0] <- "Upregulated"
df_M2A$Direction[df_M2A$padj < 0.05 & df_M2A$log2FoldChange < 0] <- "Downregulated"
df_M2A <- as.data.frame(res_Crp_pre_LFC) %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(
    Direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    Shape_Group = ifelse(abs(log2FoldChange) > window_xlimit1 | -log10(padj) > window_ylimit1, "Outlier", "Normal"),
    plot_x = case_when(
      log2FoldChange > window_xlimit1 ~ window_xlimit1,
      log2FoldChange < -window_xlimit1 ~ -window_xlimit1,
      TRUE ~ log2FoldChange
    ),
    plot_y = case_when(
      -log10(padj) > window_ylimit1 ~ window_ylimit1,
      TRUE ~ -log10(padj)),
  ) %>%
  mutate(
    plot_x = ifelse(Shape_Group == "Outlier", plot_x + runif(n(), -0.05, 0.05), plot_x),
    plot_y = ifelse(Shape_Group == "Outlier", plot_y + runif(n(), -0.05, 0.05), plot_y)
  )

M2A <- ggplot(df_M2A, aes(x = plot_x, y = plot_y, color = Direction, shape = Shape_Group)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "gray", "Upregulated" = "red")) +
  coord_cartesian(xlim = c(-window_xlimit1, window_xlimit1), ylim = c(0, window_ylimit1)) +  labs(title = "Csp2 vs. Csp5", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  annotate("text", x = -window_xlimit1 + 0.01, y = 0.2, label = "n = 11", color = "blue", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = window_xlimit1 - 0.01, y = 0.2, label = "n = 48", color = "red", hjust = 1, size = 4, fontface = "bold") +
  theme_bw() + 
  theme(legend.position = "none")
M2A

window_xlimit2 <- 1.7
window_ylimit2 <- 3
df_M2B <- as.data.frame(res_LH_pre_LFC)
df_M2B$Direction <- "Not Significant"
df_M2B$Direction[df_M2B$padj < 0.05 & df_M2B$log2FoldChange > 0] <- "Upregulated"
df_M2B$Direction[df_M2B$padj < 0.05 & df_M2B$log2FoldChange < 0] <- "Downregulated" 
df_M2B <- as.data.frame(res_LH_pre_LFC) %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(
    Direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    Shape_Group = ifelse(abs(log2FoldChange) > window_xlimit2 | -log10(padj) > window_ylimit2, "Outlier", "Normal"),
    plot_x = case_when(
      log2FoldChange > window_xlimit2 ~ window_xlimit2,
      log2FoldChange < -window_xlimit2 ~ -window_xlimit2,
      TRUE ~ log2FoldChange
    ),
    plot_y = case_when(
      -log10(padj) > window_ylimit2 ~ window_ylimit2,
      TRUE                           ~ -log10(padj)
    ),
  ) %>%
  mutate(
    plot_x = ifelse(Shape_Group == "Outlier", plot_x + runif(n(), -0.05, 0.05), plot_x),
    plot_y = ifelse(Shape_Group == "Outlier", plot_y + runif(n(), -0.05, 0.05), plot_y)
  )

M2B <- ggplot(df_M2B, aes(x = plot_x, y = -log10(padj), color = Direction, shape = Shape_Group)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "gray", "Upregulated" = "red")) +
  coord_cartesian(xlim = c(-window_xlimit2, window_xlimit2), ylim = c(0, window_ylimit2)) +  labs(title = "Csp2 vs. Csp5", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  labs(title = "Low vs. High", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  annotate("text", x = -window_xlimit2 + 0.01, y = 0.2, label = "n = 1,138", color = "blue", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = window_xlimit2 - 0.01, y = 0.2, label = "n = 420", color = "red", hjust = 1, size = 4, fontface = "bold") +
  theme_bw() + 
  theme(legend.position = "none")
M2B

#
vsd_pre <- vst(dds_pre, blind=FALSE)
ntd_pre <- normTransform(dds_pre)
select_pre <- order(rowMeans(counts(dds_pre,normalized=TRUE)),
                    decreasing=TRUE)[1:20]
df_pre <- as.data.frame(colData(dds_pre)[,c("HR_Cat")])
#pheatmap(assay(vsd_pre)[select_pre,], cluster_rows = FALSE, show_rownames = FALSE,
cluster_cols = FALSE, annotation_col = df_pre)

dev.off()
sampleDists_pre <- dist(t(assay(vsd_pre)))
sampleDistMatrix_pre <- as.matrix(sampleDists_pre)
rownames(sampleDistMatrix_pre) <- paste(vsd_pre$samples_pre.Ctrl_Exp, vsd_pre$samples_pre.SampleID, vsd_pre$samples_pre.HR_Cat, sep="-")
colnames(sampleDistMatrix_pre) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_pre,
         clustering_distance_rows=sampleDists_pre,
         clustering_distance_cols=sampleDists_pre,
         col=colors)

plotPCA(vsd_pre, intgroup=c("HR_Cat")) + theme_bw() + stat_ellipse() 
pre_PCA <- plotPCA(vsd_pre, intgroup=c("HR_Cat", "SampleID", "Cryptic", "AverageDaysinHS"), returnData = TRUE)
percentVar_pre <- round(100 * attr(pre_PCA, "percentVar"))
M2C <- ggplot(pre_PCA, aes(x=PC1, y = PC2, fill=HR_Cat, color=HR_Cat)) +
  geom_point(data=pre_PCA, size=3, aes(fill=HR_Cat, color=HR_Cat, shape=Cryptic)) +
  scale_shape_manual(values=c(17,19,13)) +
  scale_fill_manual(values = c("orange4", "dodgerblue", "goldenrod", "firebrick")) +
  scale_color_manual(values = c("orange4", "dodgerblue", "goldenrod", "firebrick")) +
  xlab(paste0("PC1: ", percentVar_pre[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_pre[2], "% variance")) +
  coord_fixed() +
  theme_bw() + stat_ellipse() + theme(axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), axis.line = element_line(colour = "black"))

#for sample ids
#ggplot(pre_PCA, aes(x=PC1, y = PC2, fill=HR_Cat, color=HR_Cat)) +
geom_point(data=pre_PCA, size=3, aes(fill=HR_Cat, color=HR_Cat)) +
  scale_fill_manual(values = c("orange4", "firebrick", "dodgerblue", "goldenrod")) +
  scale_color_manual(values = c("orange4", "firebrick", "dodgerblue", "goldenrod")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text_repel(size=3, aes(label=SampleID), max.overlaps = Inf) +
  theme_bw() + stat_ellipse() + theme(axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), legend.title = element_blank(), axis.line = element_line(colour = "black"))
#plotPCA(vsd_pre, intgroup=c("HR_Cat", "SampleID", "Cryptic", "AverageDaysinHS")) + theme_bw() 

#bar plot of DEG counts between each comparison
quartz(w=6,h=5)
M2AB <- plot_grid(M2A, M2B, ncol = 2, labels = c("A", "B"), label_size = 12, align = "h",  axis = "tb")
M2Cp <- plot_grid(NULL,M2C,NULL, ncol = 3, rel_widths = c(0.05,1,0.05))
plot_grid(M2AB,M2Cp, ncol = 1, labels = c("", "C"), label_size = 12, rel_heights = c(0.9, 1.1), align="v",axis="lr")
quartz.save("./Figures/Fig2_final.pdf", type="pdf")

####WGCNA Analysis Pre-Heat Stress Samples####
library(tseries)
wpn_pre_vsd <- getVarianceStabilizedData(dds_pre)
write.csv(wpn_pre_vsd, file = "./R_Files/vsd_pre.csv")
rv_wpn_pre <- rowVars(wpn_pre_vsd)
summary(rv_wpn_pre)
#q75_wpn_pre <- quantile( rowVars(wpn_pre_vsd), .75)
q95_wpn_pre <- quantile( rowVars(wpn_pre_vsd), .95)
expr_normalized_pre <- wpn_pre_vsd[ rv_wpn_pre > q95_wpn_pre, ]
expr_normalized_pre[1:5,1:10]
expr_normalized_pre <- expr_normalized_pre %>%
  as.data.frame() %>%
  rename("sample1" = "M_001A",
         "sample2" = "M_001B",
         "sample3" = "M_002A",
         "sample4" = "M_002B",
         "sample5" = "M_005",
         "sample6" = "H_006",
         "sample7" = "D_007A",
         "sample8" = "D_007B",
         "sample9" = "H_008A",
         "sample10" = "H_008B",
         "sample11" = "M_009",
         "sample12" = "M_010A",
         "sample13" = "M_010B",
         "sample14" = "H_011",
         "sample15" = "D_012",
         "sample16" = "M_013",
         "sample17" = "D_018",
         "sample18" = "D_020",
         "sample19" = "L_7100A",
         "sample20" = "L_7100B",
         "sample21" = "D_7105",
         "sample22" = "D_7108",
         "sample23" = "L_7112",
         "sample24" = "M_734A",
         "sample25" = "M_734B",
         "sample26" = "M_740A",
         "sample27" = "M_740B",
         "sample28" = "L_744A",
         "sample29" = "L_744B",
         "sample30" = "L_765A",
         "sample31" = "L_765B",
         "sample32" = "D_776A",
         "sample33" = "D_776B",
         "sample34" = "M_777A",
         "sample35" = "M_777B",
         "sample36" = "H_9122A",
         "sample37" = "H_9122B",
         "sample38" = "H_915A",
         "sample39" = "H_915B",
         "sample40" = "D_9162A",
         "sample41" = "D_9162B",
         "sample42" = "D_916",
         "sample43" = "H_938A",
         "sample44" = "H_938B",
         "sample45" = "L_946",
         "sample46" = "D_948A",
         "sample47" = "D_948B",
         "sample48" = "L_950",
         "sample49" = "L_952A",
         "sample50" = "L_952B",
         "sample51" = "L_956A",
         "sample52" = "L_956B",
         "sample53" = "H_960A",
         "sample54" = "H_960B",
         "sample55" = "L_994A",
         "sample56" = "L-994B",
         "sample57" = "H_999A",
         "sample58" = "H_999B") %>%
  as.matrix()

expr_normalized_pre_df <- data.frame(expr_normalized_pre) %>%
  mutate(Gene_id = row.names(expr_normalized_pre)) %>%
  pivot_longer(-Gene_id)

input_mat = t(expr_normalized_pre)
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,
  #blockSize = 30,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.85, col = "blue", lty = 2)
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#options: 9, 10, 12

picked_power = 10
temp_cor <- WGCNA::cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 8000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

table(netwk$colors) 

cor <- temp_cor
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
write_delim(module_df,
            file = "./R_Files/gene_modules_pre_new.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

dev.off()
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships: All samples pre-heat stress", y = "Modules", fill="corr") + theme(title = element_text(size=24), axis.text = element_text(size=24), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_blank(), axis.line = element_line(colour = "black"))

modules_of_interest = c("turquoise", "blue", "brown")
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id

subexpr = expr_normalized_pre[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name,y=value,group=gene_id))+
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(~module) +
  labs(title = "Module-trait Relationships: All heated samples pre-heat stress ",x = "sample ID", y = "normalized expression") + theme(title = element_text(size=24), axis.text = element_text(size=10), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_blank(), axis.line = element_line(colour = "black"))
write.csv(submod_df, file = "./R_Files/expr_normalized_pre_test.csv")


####Pre-heat Expression data based on days to bleach####
expression <- read_csv("./R_Files/expr_normalized_pre_test.csv", na = "NA") %>%
  dplyr::rename(coral_id = name)
coral_data <- read_excel("./coral_list.xlsx", na = "NA") %>%
  mutate(cryptic_species = factor(cryptic_species),
         HR_Cat = factor(HR_Cat, levels = c("high", "moderate", "low", "dead")))

gene_data <- read_excel("./coral_list.xlsx", 2, na = "NA")

expression_data <- expression %>%
  left_join(.,coral_data,by="coral_id")

#module 1 brown#
expression_data_brown <- expression_data %>%
  group_by(module) %>%
  filter(module=="brown")

M3A <- ggplot(expression_data_brown, aes(x=days_to_bleach,y=value)) +
  geom_point(color="darkorange3") +
  geom_smooth(method="lm") +
  labs(x = "", y = "Module 1\nNormalized\nExpression (vst)") +
  annotate("text", x=2, y=13, label="atop(p < 2.2*e^{-16}, R^2 == 0.025)", size=3, parse=TRUE) +
  theme_bw() + theme(legend.position="none", axis.title = element_text(size=12),axis.text = element_text(size=12)) 

module1_pre_lm <- lm(value ~ days_to_bleach, data=expression_data_brown)
summary(module1_pre_lm)
coef(module1_pre_lm)

#Exclude Module 1 genes with negative slopes to remove from the biomarkers list
#find slopes for each sample
module1slopes <- setDT(expression_data_brown)[,.(slope=coef(lm(value ~ days_to_bleach))["days_to_bleach"]),by=c("gene_id")]
#Seven genes with negative slopes
#Amillepora37638
#Amillepora17725
#Amillepora01168
#Amillepora06944
#Amillepora02439
#Amillepora32730
#Amillepora33767

#module 2 blue#
expression_data_blue <- expression_data %>%
  group_by(module) %>%
  filter(module=="blue")

M3B <- ggplot(expression_data_blue, aes(x=days_to_bleach,y=value)) +
  geom_point(color = "purple3") +
  geom_smooth(method="lm") +
  labs(x = "", y = "Module 2\nNormalized\nExpression (vst)") +
  annotate("text", x=2, y=14, label="atop(p < 2.2*e^{-16}, R^2 == 0.053)", size=3, parse=TRUE) +
  theme_bw() + theme(legend.position="none", axis.title = element_text(size=12),axis.text = element_text(size=12)) 

module2_pre_lm <- lm(value ~ days_to_bleach, data=expression_data_blue)
summary(module2_pre_lm)
coef(module2_pre_lm)

#Exclude Module 2 genes with negative slopes to remove from the biomarkers list
#find slopes for each sample
module2slopes <- setDT(expression_data_blue)[,.(slope=coef(lm(value ~ days_to_bleach))["days_to_bleach"]),by=c("gene_id")]
#Two genes with negative slopes
#Amillepora36935
#Amillepora37142

#module 3 turquoise#
expression_data_turquoise <- expression_data %>%
  group_by(module) %>%
  filter(module=="turquoise")

M3C <- ggplot(expression_data_turquoise, aes(x=days_to_bleach,y=value)) +
  geom_point(color="forestgreen") +
  geom_smooth(method="lm") +
  labs(x = "Days to Moderately Bleach", y = "Module 3\nNormalized\nExpression (vst)") +
  annotate("text", x=2, y=13.5, label="atop(p < 2.2*e^{-16}, R^2 == 0.057)", size=3, parse=TRUE) +
  theme_bw() + theme(legend.position="none", axis.title = element_text(size=12),axis.text = element_text(size=12)) 

module3_pre_lm <- lm(value ~ days_to_bleach, data=expression_data_turquoise)
summary(module3_pre_lm)
coef(module3_pre_lm)

#Exclude Module 3 genes with positive slopes to remove from the biomarkers list
#find slopes for each sample
module3slopes <- setDT(expression_data_turquoise)[,.(slope=coef(lm(value ~ days_to_bleach))["days_to_bleach"]),by=c("gene_id")]
#Five genes with positive slopes
#Amillepora20552
#Amillepora16514
#Amillepora17990
#Amillepora20781
#Amillepora33900

###Anovas for Boxplots
# === Factorial Regression: Two-Way ANOVA ===
# Check interaction between cryptic_species and HR_Cat on gene expression
cleaned_expression <- expression_data %>%
  left_join(coral_data, by = c("coral_id", "HR_Cat", "days_to_bleach", "cryptic_species"), relationship = "many-to-many") %>%
  mutate(module = recode_values(module, "brown" ~ "Module 1", "blue" ~ "Module 2", "turquoise" ~ "Module 3"),
         module = factor(module, levels = c("Module 1", "Module 2", "Module 3")),
         HR_Cat = factor(HR_Cat, labels = c("High", "Moderate", "Low", "Dead")),
         cryptic_species = factor(cryptic_species, labels = c("2", "5"))
  ) 

aov_factorial <- aov(value ~ cryptic_species * HR_Cat, data = cleaned_expression)
summary(aov_factorial)

# Post hoc test for combinations (all modules combined)
library(emmeans)
library(multcomp)
library(purrr)
emm <- emmeans(aov_factorial, ~ cryptic_species * HR_Cat)
pairs(emm, adjust = "tukey")

# List of unique modules
modules <- unique(cleaned_expression$module)

# Loop through each module and perform analysis
results_list <- lapply(modules, function(mod) {
  mod_data <- cleaned_expression %>% filter(module == mod)
  aov_mod <- aov(value ~ cryptic_species * HR_Cat, data = mod_data)
  emm_mod <- emmeans(aov_mod, ~ cryptic_species * HR_Cat)
  tukey_mod <- pairs(emm_mod, adjust = "tukey")
  cat("\n=== Module:", mod, "===\n")
  print(summary(tukey_mod))
  # NEW: Automatically generate the letters for this specific module
  cld_output <- cld(emm_mod, alpha = 0.05, Letters = letters, adjust = "tukey") %>%
    as_tibble() %>%
    mutate(
      module = mod,                  
      .group = trimws(.group),       
      y_position = max(mod_data$value, na.rm = TRUE) + 0.5 
    )
  cat("\n=== Module:", mod, "===\n")
  print(summary(tukey_mod))
  return(list(module = mod, aov = aov_mod, emmeans = emm_mod, tukey = tukey_mod, letters = cld_output))
})

all_letters <- map_dfr(results_list, ~ .x$letters)

###Boxplot including cryptic species information
#module 1 = 352, module 2 = 426, module 3 = 527

cleaned_expression_m1 <- cleaned_expression %>%
  group_by(module) %>%
  filter(module=="Module 1")

m1_letters <- results_list[[which(sapply(results_list, function(x) x$module == "Module 1"))]]$letters

M3D <- ggplot(cleaned_expression_m1, aes(x = cryptic_species, y = value, fill = HR_Cat)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
  geom_text(data = m1_letters, aes(x = cryptic_species, y = y_position, label = .group, group = HR_Cat),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 3,inherit.aes = FALSE) +
  scale_fill_manual(values = c("High" = "firebrick", "Moderate" = "goldenrod", "Low" = "dodgerblue", "Dead" = "orange4")) +
  labs(x = "", y = "", fill = "Heat Resistance") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 12),axis.title = element_text(size = 12),strip.text = element_text(size = 12),legend.text = element_text(size = 12))

cleaned_expression_m2 <- cleaned_expression %>%
  group_by(module) %>%
  filter(module=="Module 2")

m2_letters <- results_list[[which(sapply(results_list, function(x) x$module == "Module 2"))]]$letters

M3E <- ggplot(cleaned_expression_m2, aes(x = cryptic_species, y = value, fill = HR_Cat)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
  geom_text(data = m2_letters, aes(x = cryptic_species, y = y_position, label = .group, group = HR_Cat),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 3,inherit.aes = FALSE) +
  scale_fill_manual(values = c("High" = "firebrick", "Moderate" = "goldenrod", "Low" = "dodgerblue", "Dead" = "orange4")) +
  labs(x = "", y = "", fill = "Heat Resistance") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),strip.text = element_text(size = 12),legend.text = element_text(size = 12))

cleaned_expression_m3 <- cleaned_expression %>%
  group_by(module) %>%
  filter(module=="Module 3")

m3_letters <- results_list[[which(sapply(results_list, function(x) x$module == "Module 3"))]]$letters

M3F <- ggplot(cleaned_expression_m3, aes(x = cryptic_species, y = value, fill = HR_Cat)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
  geom_text(data = m3_letters, aes(x = cryptic_species, y = y_position, label = .group, group = HR_Cat),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 3,inherit.aes = FALSE) +
  scale_fill_manual(values = c("High" = "firebrick", "Moderate" = "goldenrod", "Low" = "dodgerblue", "Dead" = "orange4")) +
  labs(x = "Cryptic Species", y = "", fill = "Heat Resistance") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 12),axis.title = element_text(size = 12),strip.text = element_text(size = 12),legend.text = element_text(size = 12))

quartz(w=8,h=6)
M3ABC <- plot_grid(M3A,M3B,M3C, ncol = 1, labels = c("A", "B", "C"), label_size=12, align="hv",axis="lr")
M3DEF <- plot_grid(M3D,M3E,M3F, ncol=1, labels = c("D","E","F"), label_size = 12, align = "hv",axis = "lr")
plot_grid(M3ABC,M3DEF,ncol=2)

quartz.save("./Figures/Fig3_final.pdf", type="pdf")

####Pre-heat Bio markers and GO enrichment Summary Figure 4A-B####
install.packages("stringr")
library(stringr)

biomarkers <- read_excel("./DEG_prevspostheatstress_suppmaterials.xlsx", 5, skip = 1, na=".") %>%
  filter(Modules!="Unsorted") %>%
  filter(Biomarker!="N") %>%
  na.omit() %>%
  separate_rows(GO_Terms, sep=";") %>%
  separate(GO_Terms, into=c("GO_ID", "Ontology", "Description"),
           sep=" \\^", extra="merge") %>%
  filter(Ontology=="biological_process") %>%
  mutate(Category = case_when(str_detect(Description,"protein folding|protein refolding|chaperone|heat shock|unfolded|endoplasmic reticulum")~ "Proteostasis",
                              str_detect(Description,"ubiquitin|proteasome|protein degradation|proteolysis|autophagy|lysosome") ~ "Protein degradation",
                              str_detect(Description,"oxidative|peroxide|hydrogen peroxide|oxygen radical|reactive oxygen")~ "Oxidative stress",
                              str_detect(Description,"dna repair|recombination")~ "DNA repair",
                              str_detect(Description,"apopt|cell death")~ "Apoptosis",
                              str_detect(Description,"mitochond|electron transport|respiration|oxidative phosphorylation|atp synth|tricarboxylic acid|TCA cycle|mitochondrial membrane") ~ "Mitochondrial function",
                              str_detect(Description,"translation|ribosome")~ "Translation",
                              str_detect(Description,"signal transduction|kinase|phosphorylation|phosphatase|MAPK|calcium|second messenger|receptor") ~ "Cell signaling",
                              str_detect(Description,"ion transport|ion channel|calcium|metal ion|sodium|potassium|transporter|homeostasis") ~ "Ion transport",
                              str_detect(Description,"vesicle|transport|golgi|endocyt")~ "Membrane trafficking",
                              str_detect(Description,"cytoskeleton|actin|microtubule|tubulin|cell junction|extracellular matrix|collagen|structural") ~ "Cell structure",
                              str_detect(Description,"transcription|gene expression|rna processing|splicing|methylation")~ "Gene regulation",
                              str_detect(Description,"lipid|fatty acid|cholesterol|phospholipid|membrane lipid") ~ "Lipid metabolism",
                              str_detect(Description,"glycolysis|glucose|carbohydrate|sugar|glycogen") ~ "Carbohydrate metabolism",
                              str_detect(Description,"metabolic process|biosynthetic process|catabolic process|metabolism") ~ "General metabolism",
                              str_detect(Description,"cell cycle|mitotic|chromosome|DNA replication|replication fork|cell proliferation|growth") ~ "Cell cycle",
                              str_detect(Description,"immune|defense|inflammatory|cytokine|interferon|complement|phagocyt") ~ "Immune response",
                              str_detect(Description,"cellular stress|stress response|response to stress|response to heat|response to temperature|adaptation|homeostasis|environmental stress") ~ "Stress response",
                              TRUE ~ "Other")) %>%
  distinct(Gene_ID,Modules,Category) %>%
  dplyr::count(Modules,Category) %>%
  group_by(Modules) %>%
  mutate(Percent = 100 * n / sum(n)) %>%
  mutate(Category = factor(Category,
                           levels = c(
                             "Proteostasis",
                             "Protein degradation",
                             "Oxidative stress",
                             "DNA repair",
                             "Apoptosis",
                             "Mitochondrial function",
                             "Translation",
                             "Cell signaling",
                             "Ion transport",
                             "Membrane trafficking",
                             "Cell structure",
                             "Gene regulation",
                             "Lipid metabolism",
                             "Carbohydrate metabolism",
                             "General metabolism",
                             "Cell cycle",
                             "Immune response",
                             "Stress response",
                             "Other"
                           )))

library(viridis)

S2 <- ggviridisS2 <- ggplot(biomarkers,aes(Category,Percent,fill=Category))+
  geom_col()+
  facet_wrap(~Modules)+
  coord_flip()+
  labs(y="Number of Genes (%)", x="GO Terms (Biological Processes)")+
  scale_fill_viridis_d(option = "plasma")+
  theme_bw()+
  theme(legend.position = "none")

quartz(w=8,h=5)
plot_grid(S2)
quartz.save("./Figures/FigS2.pdf", type = "pdf")

#TopGO Pre Heat Stress Analysis#
#GO Enrichment Analysis Pre-Heat Stress
library(topGO)

#Load Trinotate annotation and prepare gene2go
anno <- read_excel("./Amillepora_trinotate_report.xlsx") %>%
  mutate(all_gos = paste(gene_ontology_blast, gene_ontology_pfam, sep = ";")) %>%
  filter(!is.na(all_gos) & all_gos != ".") %>%
  mutate(go_terms = str_extract_all(all_gos, "\\d{7}")) %>%
  filter(lengths(go_terms) > 0) %>%
  mutate(go_terms = sapply(go_terms, function(x) paste0("GO:", x, collapse = ","))) %>%
  dplyr::select(`#gene_id`, go_terms) %>%
  dplyr::rename(gene_id = `#gene_id`)

write.table(anno,
            file = "./R_files/gene2go.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

gene2GO <- readMappings(file = "./R_Files/gene2go.txt")

gene_data <- read_excel("./R_Files/gene_modules_pre_new.xlsx") %>%
  mutate(gene_id_clean = gsub("-RA$", "", gene_id)) %>%
  dplyr::select(gene_id_clean, colors) %>%
  distinct()

#Module 1 genes (Brown)
geneUniverseM1 <- gene_data$gene_id_clean[gene_data$gene_id_clean %in% names(gene2GO)]

interesting_genesM1 <- gene_data %>%
  filter(colors == "brown") %>%
  filter(gene_id_clean %in% geneUniverseM1) %>%
  pull(gene_id_clean)

geneFactorM1 <- factor(as.integer(geneUniverseM1 %in% interesting_genesM1))
names(geneFactorM1) <- geneUniverseM1

GOdata_BP_M1 <- new("topGOdata",
                    description = "GO Enrichment - BP",
                    ontology = "BP",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM1,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_BP_M1 <- runTest(GOdata_BP_M1, algorithm = "classic", statistic = "fisher")
topNodes_BP <- sum(score(result_BP_M1) < 0.05)
allRes_BP_M1 <- GenTable(GOdata_BP_M1, classicFisher = result_BP_M1,
                         orderBy = "classicFisher", topNodes = topNodes_BP)

GOdata_MF_M1 <- new("topGOdata",
                    description = "GO Enrichment - MF",
                    ontology = "MF",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM1,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_MF_M1 <- runTest(GOdata_MF_M1, algorithm = "classic", statistic = "fisher")
topNodes_MF <- sum(score(result_MF_M1) < 0.05)
allRes_MF_M1 <- GenTable(GOdata_MF_M1, classicFisher = result_MF_M1,
                         orderBy = "classicFisher", topNodes = topNodes_MF)

GOdata_CC_M1 <- new("topGOdata",
                    description = "GO Enrichment - CC",
                    ontology = "CC",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM1,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_CC_M1 <- runTest(GOdata_CC_M1, algorithm = "classic", statistic = "fisher")
topNodes_CC <- sum(score(result_CC_M1) < 0.05)
allRes_CC_M1 <- GenTable(GOdata_CC_M1, classicFisher = result_CC_M1,
                         orderBy = "classicFisher", topNodes = topNodes_CC)

#Combine and adjust p-values
allRes_BP_M1$Module <- "Module 1"; allRes_BP_M1$Ontology <- "BP"
allRes_MF_M1$Module <- "Module 1"; allRes_MF_M1$Ontology <- "MF"
allRes_CC_M1$Module <- "Module 1"; allRes_CC_M1$Ontology <- "CC"

allResM1 <- bind_rows(allRes_BP_M1, allRes_MF_M1, allRes_CC_M1)
allResM1$classicFisher <- as.numeric(allResM1$classicFisher)
allResM1$resultFisher_adj <- p.adjust(allResM1$classicFisher, method = "fdr")

write.csv(allResM1, "./R_Files/topGO_results_module1.csv", row.names = FALSE)

#Combine gene–GO pairs across all ontologies
go2genes_BP <- genesInTerm(GOdata_BP_M1)
go2genes_MF <- genesInTerm(GOdata_MF_M1)
go2genes_CC <- genesInTerm(GOdata_CC_M1)

go2genesM1 <- c(go2genes_BP, go2genes_MF, go2genes_CC)

go_gene_df_M1 <- lapply(names(go2genesM1), function(go_id) {
  genes <- go2genesM1[[go_id]]
  genes <- genes[genes %in% interesting_genesM1]
  if (length(genes) > 0) {
    data.frame(GO.ID = go_id, gene_id = genes)
  } else {
    NULL
  }
}) %>% bind_rows()

#Merge and retain all entries, no significance filtering
allResM1_with_genes_all <- go_gene_df_M1 %>%
  left_join(allResM1, by = "GO.ID") %>%
  relocate(gene_id, .after = GO.ID) %>%
  filter(!is.na(Term)) %>%
  distinct(gene_id, GO.ID, Term, Ontology, Module, .keep_all = TRUE) %>%
  arrange(classicFisher)

write.csv(allResM1_with_genes_all,
          "./R_Files/topGO_results_module1_all_entries.csv",
          row.names = FALSE)

#Module 2 genes (Blue)
geneUniverseM2 <- gene_data$gene_id_clean[gene_data$gene_id_clean %in% names(gene2GO)]

interesting_genesM2 <- gene_data %>%
  filter(colors == "blue") %>%
  filter(gene_id_clean %in% geneUniverseM2) %>%
  pull(gene_id_clean)

geneFactorM2 <- factor(as.integer(geneUniverseM2 %in% interesting_genesM2))
names(geneFactorM2) <- geneUniverseM2

GOdata_BP_M2 <- new("topGOdata",
                    description = "GO Enrichment - BP",
                    ontology = "BP",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM2,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_BP_M2 <- runTest(GOdata_BP_M2, algorithm = "classic", statistic = "fisher")
topNodes_BP <- sum(score(result_BP_M2) < 0.05)
allRes_BP_M2 <- GenTable(GOdata_BP_M2, classicFisher = result_BP_M2,
                         orderBy = "classicFisher", topNodes = topNodes_BP)

GOdata_MF_M2 <- new("topGOdata",
                    description = "GO Enrichment - MF",
                    ontology = "MF",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM2,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_MF_M2 <- runTest(GOdata_MF_M2, algorithm = "classic", statistic = "fisher")
topNodes_MF <- sum(score(result_MF_M2) < 0.05)
allRes_MF_M2 <- GenTable(GOdata_MF_M2, classicFisher = result_MF_M2,
                         orderBy = "classicFisher", topNodes = topNodes_MF)

GOdata_CC_M2 <- new("topGOdata",
                    description = "GO Enrichment - CC",
                    ontology = "CC",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM2,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_CC_M2 <- runTest(GOdata_CC_M2, algorithm = "classic", statistic = "fisher")
topNodes_CC <- sum(score(result_CC_M2) < 0.05)
allRes_CC_M2 <- GenTable(GOdata_CC_M2, classicFisher = result_CC_M2,
                         orderBy = "classicFisher", topNodes = topNodes_CC)

#Combine and adjust p-values
allRes_BP_M2$Module <- "Module 2"; allRes_BP_M2$Ontology <- "BP"
allRes_MF_M2$Module <- "Module 2"; allRes_MF_M2$Ontology <- "MF"
allRes_CC_M2$Module <- "Module 2"; allRes_CC_M2$Ontology <- "CC"

allResM2 <- bind_rows(allRes_BP_M2, allRes_MF_M2, allRes_CC_M2)
allResM2$classicFisher <- as.numeric(allResM2$classicFisher)
allResM2$resultFisher_adj <- p.adjust(allResM2$classicFisher, method = "fdr")

write.csv(allResM2, "./R_Files/topGO_results_module2.csv", row.names = FALSE)

#Combine gene–GO pairs across all ontologies
go2genes_BP <- genesInTerm(GOdata_BP_M2)
go2genes_MF <- genesInTerm(GOdata_MF_M2)
go2genes_CC <- genesInTerm(GOdata_CC_M2)

go2genesM2 <- c(go2genes_BP, go2genes_MF, go2genes_CC)

go_gene_df_M2 <- lapply(names(go2genesM2), function(go_id) {
  genes <- go2genesM2[[go_id]]
  genes <- genes[genes %in% interesting_genesM2]
  if (length(genes) > 0) {
    data.frame(GO.ID = go_id, gene_id = genes)
  } else {
    NULL
  }
}) %>% bind_rows()

#Merge and retain all entries, no significance filtering
allResM2_with_genes_all <- go_gene_df_M2 %>%
  left_join(allResM2, by = "GO.ID") %>%
  relocate(gene_id, .after = GO.ID) %>%
  filter(!is.na(Term)) %>%
  distinct(gene_id, GO.ID, Term, Ontology, Module, .keep_all = TRUE) %>%
  arrange(classicFisher)

write.csv(allResM2_with_genes_all,
          "./R_Files/topGO_results_module2_all_entries.csv",
          row.names = FALSE)

#Module 1 genes (Turquoise)
geneUniverseM3 <- gene_data$gene_id_clean[gene_data$gene_id_clean %in% names(gene2GO)]

interesting_genesM3 <- gene_data %>%
  filter(colors == "turquoise") %>%
  filter(gene_id_clean %in% geneUniverseM3) %>%
  pull(gene_id_clean)

geneFactorM3 <- factor(as.integer(geneUniverseM3 %in% interesting_genesM3))
names(geneFactorM3) <- geneUniverseM3

GOdata_BP_M3 <- new("topGOdata",
                    description = "GO Enrichment - BP",
                    ontology = "BP",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM3,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_BP_M3 <- runTest(GOdata_BP_M3, algorithm = "classic", statistic = "fisher")
topNodes_BP <- sum(score(result_BP_M3) < 0.05)
allRes_BP_M3 <- GenTable(GOdata_BP_M3, classicFisher = result_BP_M3,
                         orderBy = "classicFisher", topNodes = topNodes_BP)

GOdata_MF_M3 <- new("topGOdata",
                    description = "GO Enrichment - MF",
                    ontology = "MF",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM3,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_MF_M3 <- runTest(GOdata_MF_M3, algorithm = "classic", statistic = "fisher")
topNodes_MF <- sum(score(result_MF_M3) < 0.05)
allRes_MF_M3 <- GenTable(GOdata_MF_M3, classicFisher = result_MF_M3,
                         orderBy = "classicFisher", topNodes = topNodes_MF)

GOdata_CC_M3 <- new("topGOdata",
                    description = "GO Enrichment - CC",
                    ontology = "CC",
                    geneSel = function(x) x == 1,
                    allGenes = geneFactorM3,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
result_CC_M3 <- runTest(GOdata_CC_M3, algorithm = "classic", statistic = "fisher")
topNodes_CC <- sum(score(result_CC_M3) < 0.05)
allRes_CC_M3 <- GenTable(GOdata_CC_M3, classicFisher = result_CC_M3,
                         orderBy = "classicFisher", topNodes = topNodes_CC)

#Combine and adjust p-values
allRes_BP_M3$Module <- "Module 3"; allRes_BP_M3$Ontology <- "BP"
allRes_MF_M3$Module <- "Module 3"; allRes_MF_M3$Ontology <- "MF"
allRes_CC_M3$Module <- "Module 3"; allRes_CC_M3$Ontology <- "CC"

allResM3 <- bind_rows(allRes_BP_M3, allRes_MF_M3, allRes_CC_M3)
allResM3$classicFisher <- as.numeric(allResM3$classicFisher)
allResM3$resultFisher_adj <- p.adjust(allResM3$classicFisher, method = "fdr")

write.csv(allResM3, "./R_Files/topGO_results_module3.csv", row.names = FALSE)

#Combine gene–GO pairs across all ontologies
go2genes_BP <- genesInTerm(GOdata_BP_M3)
go2genes_MF <- genesInTerm(GOdata_MF_M3)
go2genes_CC <- genesInTerm(GOdata_CC_M3)

go2genesM3 <- c(go2genes_BP, go2genes_MF, go2genes_CC)

go_gene_df_M3 <- lapply(names(go2genesM3), function(go_id) {
  genes <- go2genesM3[[go_id]]
  genes <- genes[genes %in% interesting_genesM3]
  if (length(genes) > 0) {
    data.frame(GO.ID = go_id, gene_id = genes)
  } else {
    NULL
  }
}) %>% bind_rows()

#Merge and retain all entries, no significance filtering
allResM3_with_genes_all <- go_gene_df_M3 %>%
  left_join(allResM3, by = "GO.ID") %>%
  relocate(gene_id, .after = GO.ID) %>%
  filter(!is.na(Term)) %>%
  distinct(gene_id, GO.ID, Term, Ontology, Module, .keep_all = TRUE) %>%
  arrange(classicFisher)

write.csv(allResM3_with_genes_all,
          "./R_Files/topGO_results_module3_all_entries.csv",
          row.names = FALSE)

#MODULE 1-3 (COMPILED)

combined_results_all <- bind_rows(
  allResM1_with_genes_all,
  allResM2_with_genes_all,
  allResM3_with_genes_all) %>%
  distinct(Module, Ontology, GO.ID, gene_id, .keep_all = TRUE) %>%
  relocate(gene_id, .after = GO.ID) %>%
  arrange(Module, Ontology, classicFisher)

write.csv(combined_results_all,"./R_Files/topGO_combinedresults.csv",row.names = FALSE)

combined_results_all %>%
  dplyr::count(Module, GO.ID, gene_id) %>%
  filter(n > 1) %>%
  dplyr::select(gene_id, Module) %>%
  distinct()

library(tidytext)

#Filter BP-only GO terms and rename modules numerically
dot_data_bp <- combined_results_all %>%
  filter(Ontology == "BP") %>%
  group_by(Module, GO.ID, Term, Ontology) %>%
  summarise(pval = min(classicFisher),
            gene_count = n_distinct(gene_id),
            .groups = "drop") %>%
  filter(pval <= 0.05)

#Select top 15 terms per module
top_15_per_module <- dot_data_bp %>%
  group_by(Module) %>%
  slice_min(pval, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Module, pval) %>%
  mutate(y = factor(Term, levels = rev(unique(Term))))

box_bounds <- top_15_per_module %>%
  mutate(y_num = as.numeric(y)) %>%
  group_by(Module) %>%
  summarise(
    ymin = min(y_num) - 0.2,
    ymax = max(y_num) + 0.2,
    .groups = "drop")

S3 <- ggplot(top_15_per_module, aes(x = 1, y = y)) +
  geom_point(aes(size = gene_count, color = pval)) +
  geom_rect(data = box_bounds,
            aes(xmin = 0.7, xmax = 1.3, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 1) +
  scale_color_viridis_c(direction = -1) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  annotate("text", x=0.8,y=44,label="Module 1")+
  annotate("text", x=0.8,y=29,label="Module 2")+
  annotate("text", x=0.8,y=14,label="Module 3")+
  labs(title = "Top 15 Biological\nProcesses GO Terms", x = NULL, y = "GO Term",color = "p-value", size = "Gene Count") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"))

quartz(w=8,h=6)
plot_grid(S3)
quartz.save("./Figures/FigS3.pdf", type = "pdf")

####POST:heated vs control####
samples_post <- read_excel(file.path(dir, "sampleIDs_data.xlsx"), 2)
files_post <- file.path(dir, "Samples_Current", samples_post$SampleID, "abundance.h5")
names(files_post) <- paste0("sample", 1:45)
all(file.exists(files_post))
txi_post <- tximport(files_post, type = "kallisto", txOut = TRUE)
head(txi_post$counts)

sampleTable_post <- data.frame(samples_post$SampleID, samples_post$HR_Cat, samples_post$Treatment, samples_post$AverageDaysinHS)
rownames(sampleTable_post) <- colnames(txi_post$counts)
all(rownames(sampleTable_post) %in% colnames(txi_post$counts))
all(rownames(sampleTable_post) == colnames(txi_post$counts))
colnames(sampleTable_post) <- c("SampleID", "HR_Cat", "Treatment", "AverageDaysinHS")
sampleTable_post$Group <- factor(paste(sampleTable_post$Treatment, sampleTable_post$HR_Cat, sep="_"))

dds_post <- DESeqDataSetFromTximport(txi_post, colData = sampleTable_post, design = ~ Group)
dds_post <- DESeq(dds_post)

resultsNames(dds_post)

#ALL heated vs control
res_ALL_H_vs_C <- results(dds_post, contrast=list(c("Group_Heated_High_vs_Control_High", "Group_Heated_Low_vs_Control_High", "Group_Heated_Moderate_vs_Control_High"), c("Group_Control_Low_vs_Control_High", "Group_Control_Moderate_vs_Control_High")),
                          listValues=c(1/3, -1/3), 
                          alpha=0.05, parallel=TRUE)
res_ALL_H_vs_C
res_ALL_H_vs_C <- res_ALL_H_vs_C[order(res_ALL_H_vs_C$pvalue),]
write.csv(res_ALL_H_vs_C,"./R_Files/ALL_treatments_post_raw.csv")
summary(res_ALL_H_vs_C)
metadata(res_ALL_H_vs_C)$filterThreshold
sum(res_ALL_H_vs_C$padj < 0.05, na.rm = TRUE)
#353 out of 27460
#LFC > 0 63
#LFC < 0 290

res_ALL_H_vs_C_LFC <- lfcShrink(dds_post, 
                                contrast=list(
                                  c("Group_Heated_High_vs_Control_High", 
                                    "Group_Heated_Low_vs_Control_High", 
                                    "Group_Heated_Moderate_vs_Control_High"), 
                                  c("Group_Control_Low_vs_Control_High", 
                                    "Group_Control_Moderate_vs_Control_High")),
                                res=res_ALL_H_vs_C, # This passes the pre-calculated listValues
                                type="ashr")
df_ALL_H_vs_C_LFC <- as.data.frame(res_ALL_H_vs_C_LFC)
df_ALL_H_vs_C_LFC <- df_ALL_H_vs_C_LFC[order(df_ALL_H_vs_C_LFC$pvalue), ]
summary(res_ALL_H_vs_C_LFC, alpha = 0.05)
sum(df_ALL_H_vs_C_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_ALL_H_vs_C_LFC, "./R_Files/ALL_treatments_post_ADJUSTED.csv")
#353 out of 27460
#LFC > 0 63
#LFC < 0 290

#Low vs Moderate Resistance (Heated Only)
res_post_LM <- results(dds_post, contrast=c("Group", "Heated_Low", "Heated_Moderate"), alpha=0.05, parallel=TRUE)
res_post_LM
res_post_LM <- res_post_LM[order(res_post_LM$pvalue),]
write.csv(res_post_LM,"./R_Files/heatedLM_post_raw.csv")
summary(res_post_LM)
metadata(res_post_LM)$filterThreshold
sum(res_post_LM$padj < 0.05, na.rm = TRUE)
#47 out of 27460
#LFC > 0 34
#LFC < 0 13

res_post_LM_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_Low", "Heated_Moderate"), res=res_post_LM, type="ashr")
df_post_LM_LFC <- as.data.frame(res_post_LM_LFC)
df_post_LM_LFC <- df_post_LM_LFC[order(df_post_LM_LFC$pvalue), ]
summary(res_post_LM_LFC, alpha = 0.05)
sum(df_post_LM_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_LM_LFC, "./R_Files/heatedLM_post_ADJUSTED.csv")
#47 out of 27460
#LFC > 0 34
#LFC < 0 13

#Low vs High Resistance (Heated Only)
res_post_LH <- results(dds_post, contrast=c("Group", "Heated_Low", "Heated_High"), alpha=0.05, parallel=TRUE)
res_post_LH
res_post_LH <- res_post_LH[order(res_post_LH$pvalue),]
write.csv(res_post_LH,"./R_Files/heatedLH_post_raw.csv")
summary(res_post_LH)
metadata(res_post_LH)$filterThreshold
sum(res_post_LH$padj < 0.05, na.rm = TRUE)
#9 out of 27460
#LFC > 0 2
#LFC < 0 7

res_post_LH_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_Low", "Heated_High"), res=res_post_LH, type="ashr")
df_post_LH_LFC <- as.data.frame(res_post_LH_LFC)
df_post_LH_LFC <- df_post_LH_LFC[order(df_post_LH_LFC$pvalue), ]
summary(res_post_LH_LFC, alpha = 0.05)
sum(df_post_LH_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_LH_LFC, "./R_Files/heatedLH_post_ADJUSTED.csv")
#9 out of 27460
#LFC > 0 2
#LFC < 0 7

#Moderate vs High Resistance (heated)
res_post_MH <- results(dds_post, contrast=c("Group", "Heated_Moderate", "Heated_High"), alpha=0.05, parallel=TRUE)
res_post_MH
res_post_MH <- res_post_MH[order(res_post_MH$pvalue),]
write.csv(res_post_MH,"./R_Files/heatedMH_post_raw.csv")
summary(res_post_MH)
metadata(res_post_MH)$filterThreshold
sum(res_post_MH$padj < 0.05, na.rm = TRUE)
#26 out of 27460
#LFC > 0 3
#LFC < 0 23

res_post_MH_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_Moderate", "Heated_High"), res=res_post_MH, type="ashr")
df_post_MH_LFC <- as.data.frame(res_post_MH_LFC)
df_post_MH_LFC <- df_post_MH_LFC[order(df_post_MH_LFC$pvalue), ]
summary(res_post_MH_LFC, alpha = 0.05)
sum(df_post_MH_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_MH_LFC, "./R_Files/heatedMH_post_ADJUSTED.csv")
#26 out of 27460
#LFC > 0 3
#LFC < 0 23

#Low Resistance (heated vs control)
res_post_HC_Low <- results(dds_post, contrast=c("Group", "Heated_Low", "Control_Low"), alpha=0.05, parallel=TRUE)
res_post_HC_Low
res_post_HC_Low <- res_post_HC_Low[order(res_post_HC_Low$pvalue),]
write.csv(res_post_HC_Low,"./R_Files/Low_HC_post_raw.csv")
summary(res_post_HC_Low)
metadata(res_post_HC_Low)$filterThreshold
sum(res_post_HC_Low$padj < 0.05, na.rm = TRUE)
#26 out of 27460
#LFC > 0 6
#LFC < 0 20

res_post_HC_Low_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_Low", "Control_Low"), res=res_post_HC_Low, type="ashr")
df_post_HC_Low_LFC <- as.data.frame(res_post_HC_Low_LFC)
df_post_HC_Low_LFC <- df_post_HC_Low_LFC[order(df_post_HC_Low_LFC$pvalue), ]
summary(res_post_HC_Low_LFC, alpha = 0.05)
sum(df_post_HC_Low_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_HC_Low_LFC, "./R_Files/Low_HC_post_ADJUSTED.csv")
#26 out of 27460
#LFC > 0 6
#LFC < 0 20

#Moderate Resistance (heated vs control)
res_post_HC_Mod <- results(dds_post, contrast=c("Group", "Heated_Moderate", "Control_Moderate"), alpha=0.05, parallel=TRUE)
res_post_HC_Mod
res_post_HC_Mod <- res_post_HC_Mod[order(res_post_HC_Mod$pvalue),]
write.csv(res_post_HC_Mod,"./R_Files/Mod_HC_post_raw.csv")
summary(res_post_HC_Mod)
metadata(res_post_HC_Mod)$filterThreshold
sum(res_post_HC_Mod$padj < 0.05, na.rm = TRUE)
#308 out of 27460
#LFC > 0 88
#LFC < 0 220

res_post_HC_Mod_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_Moderate", "Control_Moderate"), res=res_post_HC_Mod, type="ashr")
df_post_HC_Mod_LFC <- as.data.frame(res_post_HC_Mod_LFC)
df_post_HC_Mod_LFC <- df_post_HC_Mod_LFC[order(df_post_HC_Mod_LFC$pvalue), ]
summary(res_post_HC_Mod_LFC, alpha = 0.05)
sum(df_post_HC_Mod_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_HC_Mod_LFC, "./R_Files/Mod_HC_post_ADJUSTED.csv")
#308 out of 27460
#LFC > 0 88
#LFC < 0 220

#High Resistance (heated vs control)
res_post_HC_High <- results(dds_post, contrast=c("Group", "Heated_High", "Control_High"), alpha=0.05, parallel=TRUE)
res_post_HC_High
res_post_HC_High <- res_post_HC_High[order(res_post_HC_High$pvalue),]
write.csv(res_post_HC_High,"./R_Files/High_HC_post_raw.csv")
summary(res_post_HC_High)
metadata(res_post_HC_High)$filterThreshold
sum(res_post_HC_High$padj < 0.05, na.rm = TRUE)
#52 out of 27460
#LFC > 0 10
#LFC < 0 42

res_post_HC_High_LFC <- lfcShrink(dds_post, contrast=c("Group", "Heated_High", "Control_High"), res=res_post_HC_High, type="ashr")
df_post_HC_High_LFC <- as.data.frame(res_post_HC_High_LFC)
df_post_HC_High_LFC <- df_post_HC_High_LFC[order(df_post_HC_High_LFC$pvalue), ]
summary(res_post_HC_High_LFC, alpha = 0.05)
sum(df_post_HC_High_LFC$padj < 0.05, na.rm = TRUE)
write.csv(df_post_HC_High_LFC, "./R_Files/High_HC_post_ADJUSTED.csv")
#52 out of 27460
#LFC > 0 10
#LFC < 0 42

####Fig. S4 MA plots Post All####
pdf("./Figures/S4_PostStress_Expression_Comparisons.pdf", width = 11, height = 8.5)
par(mfrow=c(2,2), mar=c(4,4,2,1), oma = c(0, 0, 3, 0))
xlim <- c(1,1e4); ylim <- c(-3,3); 
plotMA(res_ALL_H_vs_C_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "All Heated vs. All Control")
plotMA(res_post_LM_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Low vs. Moderate")
plotMA(res_post_LH_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Low vs. High")
plotMA(res_post_MH_LFC, xlim=xlim, ylim=ylim, alpha=0.05, cex=1.5, main = "Moderate vs. High")

mtext("Post-Heat Stress Comparisons Between Heat Resistance Categories", 
      outer = TRUE, 
      side = 3, 
      cex = 1.5, 
      font = 2, 
      line = 1)

dev.off()

####Fig. 4A&B Post-Stress PCA and Differential Expression Counts####
vsd_post <- vst(dds_post, blind=FALSE)
ntd_post <- normTransform(dds_post)
select_post <- order(rowMeans(counts(dds_post,normalized=TRUE)),
                     decreasing=TRUE)[1:45]
df_post <- as.data.frame(colData(dds_post)[,c("HR_Cat")])
pheatmap(assay(vsd_post)[select_post,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df_post)

dev.off()
sampleDists_post <- dist(t(assay(vsd_post)))
sampleDistMatrix_post <- as.matrix(sampleDists_post)
rownames(sampleDistMatrix_post) <- paste(vsd_post$samples_post.Treatment, vsd_post$samples_post.SampleID, vsd_post$samples_post.HR_Cat, sep="-")
colnames(sampleDistMatrix_post) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_post,
         clustering_distance_rows=sampleDists_post,
         clustering_distance_cols=sampleDists_post,
         col=colors)

plotPCA(vsd_post, intgroup=c("HR_Cat")) + theme_bw() + stat_ellipse() 
post_PCA <- plotPCA(vsd_post, intgroup=c("HR_Cat", "SampleID", "Treatment"), returnData = TRUE)
percentVar_post <- round(100 * attr(post_PCA, "percentVar"))
M4A <- ggplot(post_PCA, aes(x=PC1, y = PC2, fill=Treatment, color=Treatment)) +
  geom_point(data=post_PCA, size=5, aes(fill=Treatment, color=Treatment)) +
  scale_shape_manual(values=c(17,19)) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) +
  scale_color_manual(values = c("dodgerblue", "firebrick")) +
  xlab(paste0("PC1: ", percentVar_post[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_post[2], "% variance")) +
  #coord_fixed() +
  theme_bw() + stat_ellipse() + theme(axis.text = element_text(size=12), axis.title = element_text(size=12), 
                                      legend.text = element_text(size=12), legend.title = element_blank(), 
                                      legend.position = c(0.85, 0.15),
                                      legend.background = element_blank(),
                                      legend.key = element_blank(),
                                      axis.line = element_line(colour = "black"))

#bar plot of DEG counts between each comparison
get_counts <- function(res_obj, comparison_name) {
  as.data.frame(res_obj) %>%
    filter(padj < 0.05) %>%
    summarise(
      Upregulated = sum(log2FoldChange > 0, na.rm = TRUE),
      Downregulated = sum(log2FoldChange < 0, na.rm = TRUE)
    ) %>%
    mutate(Comparison = comparison_name)
}

plot_data <- plot_data %>%
  mutate(Comparison = case_when(
    Comparison == "All Heated vs. Control" ~ "All Heated / Control",
    Comparison == "Low vs. Moderate"       ~ "Low / Moderate",
    Comparison == "Low vs. High"           ~ "Low / High",
    Comparison == "Moderate vs. High"      ~ "Moderate / High"
  ))

plot_data$Comparison <- factor(plot_data$Comparison, 
                               levels = c("All Heated / Control", 
                                          "Low / Moderate", 
                                          "Low / High", 
                                          "Moderate / High"))

M4B <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black", width = 0.7) +
  geom_text(aes(label = Count), position = position_dodge(0.7), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Upregulated" = "#de2d26", "Downregulated" = "#3182bd")) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  annotate("text", x = 0.5, y = max(plot_data$Count) * 1.30, label = "Summary of DEG Counts", fontface = "bold", hjust = 0, size = 4) +
  annotate("text", x = 0.5, y = max(plot_data$Count) * 1.20, label = "Direction: Group 1 / Group 2", fontface = "italic", hjust = 0, size = 3, color = "gray30") +
  labs(x = NULL, y = "Number of Genes", fill = NULL) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = c(0.75, 0.60), legend.background = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

#Summary of go terms post-heat stress all heated vs all control
deg <- read_excel("./R_Files/ALL_treatments_post_ADJUSTED.xlsx", 2, na=".") %>%
  na.omit() %>%
  separate_rows(goterm, sep=";") %>%
  separate(goterm,
           into=c("GO_ID","Ontology","Description"),
           sep=" \\^",
           extra="merge") %>%
  filter(Ontology == "biological_process") %>%
  mutate(Description = tolower(Description)) %>%
  mutate(Category = case_when(
    str_detect(Description, "negative regulation|decrease|down-regulation|inhibition") ~ "Other",
    str_detect(Description, "dna repair|dna damage|double-strand|recombination|chromatin|histone|dna replication|dna methyl") ~ "DNA repair & genome maintenance", 
    str_detect(Description, "protein folding|unfolded|endoplasmic reticulum|erad|ubiquitin|proteasome") ~ "Proteostasis", 
    str_detect(Description, "oxidative|hydrogen peroxide|peroxide|catalase|oxylipin") ~ "Oxidative stress", 
    str_detect(Description, "lipid|fatty acid|sterol|cholesterol|phospholipid|membrane repair|golgi|vesicle") ~ "Lipid & membrane homeostasis", 
    str_detect(Description, "translation|ribosome|transcription|rna|splicing") ~ "Translation & RNA metabolism", 
    str_detect(Description, "signal|mapk|calcium|smad|fgf|kinase") ~ "Cell signaling", 
    TRUE ~ "Other")) %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>% 
  distinct(AmilGene, Direction, Category) %>% 
  dplyr::count(Direction, Category) %>% 
  filter(Category != "Other") %>% 
  group_by(Direction) %>% 
  mutate(
    Proportion = n / sum(n),
    Percentage = round(Proportion * 100, 1)
  ) %>% 
  ungroup()

plot_data <- data.frame(
  Direction = c(rep("Downregulated", 6),rep("Upregulated", 6)),
  Category = c(
    "Cell signaling", "DNA repair & genome maintenance", "Lipid & membrane homeostasis", 
    "Oxidative stress", "Proteostasis", "Translation & RNA metabolism",
    "Cell signaling", "DNA repair & genome maintenance", "Lipid & membrane homeostasis", 
    "Oxidative stress", "Proteostasis", "Translation & RNA metabolism"),
  n = c(34, 9, 16, 6, 15, 35,  10, 8, 8, 2, 2, 10),
  Percentage = c(29.6, 7.8, 13.9, 5.2, 13.0, 30.4,25.0, 20.0, 20.0, 5.0, 5.0, 25.0))

plot_data <- plot_data %>%
  mutate(Category = factor(Category, levels = c(
    "Translation & RNA metabolism",
    "Cell signaling",
    "Lipid & membrane homeostasis",
    "Proteostasis",
    "DNA repair & genome maintenance",
    "Oxidative stress"
  )))

M4C <- ggplot(plot_data, aes(x = Category, y = Percentage, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_text(
    aes(label = paste0(Percentage, "% (n=", n, ")")), 
    position = position_dodge(width = 0.8), 
    hjust = -0.1, 
    size = 2.8, 
    fontface = "plain"
  ) +
  scale_fill_manual(values = c("Downregulated" = "#4682B4", "Upregulated" = "#CD5C5C")) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  coord_flip() + 
  theme_bw() + 
  labs(
    x = "Functional Category",
    y = "Functional Terms (%)",
    fill = NULL) +
  theme(
    axis.text.x = element_text(color = "black", face = "plain"), 
    axis.text.y = element_text(color = "black", face = "plain"),
    axis.ticks = element_line(color = "black"),
    legend.position = "top")

quartz(w=8.5,h=6)
M4AB <- plot_grid(M4A,M4B, labels = c("A","B"), label_size = 12, ncol = 2, rel_widths = c(1.7, 1),align = "h",axis = "bt")
plot_grid(M4AB,M4C,labels = c("", "C"),label_size = 12, ncol = 1, rel_heights = c(1,1.2))
quartz.save("./Figures/Fig4_final.pdf", type = "pdf")

####WGCNA Analysis Post-Heat Stress Samples####
library(tseries)
wpn_post_vsd <- getVarianceStabilizedData(dds_post)
write.csv(wpn_post_vsd, file = "./R_Files/vsd_post.csv")
rv_wpn_post <- rowVars(wpn_post_vsd)
summary(rv_wpn_post)
q75_wpn_post <- quantile( rowVars(wpn_post_vsd), .75)
q95_wpn_post <- quantile( rowVars(wpn_post_vsd), .95)
expr_normalized_post <- wpn_post_vsd[ rv_wpn_post > q95_wpn_post, ]
expr_normalized_post[1:5,1:10]
expr_normalized_post <- expr_normalized_post %>%
  as.data.frame() %>%
  rename("sample1" = "HM_001",
         "sample2" = "CM_002",
         "sample3" = "HM_002",
         "sample4" = "CH_003",
         "sample5" = "HH_003",
         "sample6" = "CM_005",
         "sample7" = "HM_005",
         "sample8" = "CH_006",
         "sample9" = "HH_006",
         "sample10" = "CH_008",
         "sample11" = "HH_008",
         "sample12" = "CM_009",
         "sample13" = "CM_010",
         "sample14" = "HM_010",
         "sample15" = "CH_011",
         "sample16" = "HH_011",
         "sample17" = "CM_013",
         "sample18" = "CL_7100",
         "sample19" = "HL_7100",
         "sample20" = "CL_7112",
         "sample21" = "HL_7112",
         "sample22" = "CM_734",
         "sample23" = "CM_740",
         "sample24" = "CL_744",
         "sample25" = "HL_744",
         "sample26" = "CL_765",
         "sample27" = "HL_765",
         "sample28" = "CM_777",
         "sample29" = "HM_777",
         "sample30" = "CH_9122",
         "sample31" = "HH_9122",
         "sample32" = "CH_915",
         "sample33" = "HH_915",
         "sample34" = "CH_938",
         "sample35" = "HH_938",
         "sample36" = "CL_946",
         "sample37" = "CL_950",
         "sample38" = "HL_950",
         "sample39" = "CL_952",
         "sample40" = "HL_952",
         "sample41" = "CL_956",
         "sample42" = "HL_956",
         "sample43" = "HH_960",
         "sample44" = "CL_994",
         "sample45" = "HL_994") %>%
  as.matrix()


expr_normalized_post_df <- data.frame(expr_normalized_post) %>%
  mutate(Gene_id = row.names(expr_normalized_post)) %>%
  pivot_longer(-Gene_id)

input_mat = t(expr_normalized_post)
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,
  #blockSize = 30,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.85, col = "blue", lty = 2)
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#options: 8, 9, 10, 12

picked_power = 9
temp_cor <- WGCNA::cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 2000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

table(netwk$colors) 

cor <- temp_cor
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
write_delim(module_df,
            file = "./R_Files/gene_modules_post.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

dev.off()
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships: All samples post-heat stress", y = "Modules", fill="corr") + theme(title = element_text(size=24), axis.text = element_text(size=24), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_blank(), axis.line = element_line(colour = "black"))

modules_of_interest = c("turquoise", "blue")
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id

subexpr = expr_normalized_post[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name,y=value,group=gene_id))+
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(~module) +
  labs(title = "Module-trait Relationships: All samples post-heat stress ",x = "sample ID", y = "normalized expression") + theme(title = element_text(size=24), axis.text = element_text(size=10), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_blank(), axis.line = element_line(colour = "black"))
write.csv(submod_df, file = "./R_Files/expr_normalized_post.csv")

####post-heat stress expression of pre-heat stress modules
expression_post <- read_csv("./R_Files/vsd_post.csv", na = "NA") %>%
  as.data.frame() %>%
  rename("sample1" = "HM_001",
         "sample2" = "CM_002",
         "sample3" = "HM_002",
         "sample4" = "CH_003",
         "sample5" = "HH_003",
         "sample6" = "CM_005",
         "sample7" = "HM_005",
         "sample8" = "CH_006",
         "sample9" = "HH_006",
         "sample10" = "CH_008",
         "sample11" = "HH_008",
         "sample12" = "CM_009",
         "sample13" = "CM_010",
         "sample14" = "HM_010",
         "sample15" = "CH_011",
         "sample16" = "HH_011",
         "sample17" = "CM_013",
         "sample18" = "CL_7100",
         "sample19" = "HL_7100",
         "sample20" = "CL_7112",
         "sample21" = "HL_7112",
         "sample22" = "CM_734",
         "sample23" = "CM_740",
         "sample24" = "CL_744",
         "sample25" = "HL_744",
         "sample26" = "CL_765",
         "sample27" = "HL_765",
         "sample28" = "CM_777",
         "sample29" = "HM_777",
         "sample30" = "CH_9122",
         "sample31" = "HH_9122",
         "sample32" = "CH_915",
         "sample33" = "HH_915",
         "sample34" = "CH_938",
         "sample35" = "HH_938",
         "sample36" = "CL_946",
         "sample37" = "CL_950",
         "sample38" = "HL_950",
         "sample39" = "CL_952",
         "sample40" = "HL_952",
         "sample41" = "CL_956",
         "sample42" = "HL_956",
         "sample43" = "HH_960",
         "sample44" = "CL_994",
         "sample45" = "HL_994")

coral_post_data <- read_excel("./coral_list.xlsx", 3, na = "NA") %>%
  mutate(coral_id=as.factor(coral_id))

expression_post_data <- expression_post %>%
  gather(coral_id, value, 2:46) %>%
  mutate(coral_id=as.factor(coral_id)) %>%
  left_join(.,coral_post_data,by="coral_id") 

coexpr_post <- read_csv("./R_Files/expr_normalized_post.csv") %>%
  dplyr::rename("coral_id"="name") %>%
  left_join(.,coral_post_data,by="coral_id") %>%
  filter(treatment=="heated") %>%
  mutate(module = case_match(module, "blue" ~"Module 4", "turquoise" ~ "Module 5"))

quartz()
S5A <- ggplot(coexpr_post, aes(x=days_to_bleach,y=value, color=module)) +
  geom_point(alpha=0.3) +
  geom_smooth(method = "lm", fill = NA) +
  scale_color_manual(
    name = "Post-Heat Stress", 
    values = c("Module 4" = "darkcyan", "Module 5" = "magenta3"),
    labels = c("Module 4" = "Module 4 (**)", "Module 5" = "Module 5 (*)")) +
  labs(x = "Heat Stress Exposure Days", y = "Post-Stress Normalized\nExpression (vst)") + 
  theme_bw() + 
  ylim(5.8,20)+
  theme(
    legend.position = c(0.7, 0.85), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12))


#post-heat stress module 4#
coexpr_post_4 <- coexpr_post %>%
  group_by(module) %>%
  filter(module=="Module 4")

module4_post_lm <- lm(value ~ days_to_bleach, data=coexpr_post_4)
summary(module4_post_lm)
#p=1.22e-10

coexpr_post_5 <- coexpr_post %>%
  group_by(module) %>%
  filter(module=="Module 5")

module5_post_lm <- lm(value ~ days_to_bleach, data=coexpr_post_5)
summary(module5_post_lm)                                                                                                                                                                                                                                                                       
#p=0.0182

#Compare post-heat stress expression against the pre-heat stress co-expression modules
pre_modules <- read_tsv("./R_Files/gene_modules_pre_new.txt") 

expression_prepost <- expression_post_data %>%
  dplyr::rename("gene_id" = 1) %>%
  left_join(.,pre_modules, by = "gene_id") %>%
  na.omit() %>%
  mutate(colors = case_match(colors, "brown" ~ "1",
                             "blue" ~ "2",
                             "turquoise" ~ "3")) %>%
  filter(treatment=="heated") %>%
  na.omit()

S5B <- ggplot(expression_prepost, aes(x = days_to_bleach, y = value, color = colors)) +
  geom_point(alpha=0.3) +
  geom_smooth(method = "lm", fill = NA) +
  scale_color_manual(
    name = "Pre-Heat Stress", 
    values = c("1" = "darkorange3", "2" = "purple3", "3" = "forestgreen"),
    labels = c("1" = "Module 1 (**)", "2" = "Module 2 (ns)", "3" = "Module 3 (**)")
  )+
  labs(
    x = "Heat Stress Exposure Days", 
    y = "") + 
  theme_bw() +
  ylim(5.8,20)+
  theme(
    legend.position = c(0.6, 0.85), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12))

separate_models <- expression_prepost %>% 
  mutate(colors = as.factor(colors)) %>% 
  group_split(colors)

walk(separate_models, function(sub_df) {
  module_number <- unique(sub_df$colors)
  cat("\n===================================\n")
  cat("SEPARATE LM FOR MODULE:", module_number, "\n")
  cat("===================================\n")
  model <- lm(value ~ days_to_bleach, data = sub_df)
  print(summary(model))
})

post_goterms <- read_excel("./DEG_prevspostheatstress_suppmaterials.xlsx", 9, na=".") %>%
  na.omit() %>%
  filter(Modules!="Unsorted")

quartz(w=8.5,h=6)
plot_grid(S5A,S5B, labels = c("A","B"))
quartz.save("./Figures/FigS5.pdf", type = "pdf")