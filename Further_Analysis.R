 

library("DESeq2") 

library("data.table") 

library(purrr) 

library(scater) 

library(cowplot) 

library(png) 

library(stringr) 

library("tidyverse") 

library("ggplot2") 

library("patchwork") 

library("dplyr") 

library("DT") 

library("Seurat") 

 

#Read in RDS object created by CTC Pipeline 

SO <- readRDS("C:/Users/changmat/Documents/CTC_Analysis_Results/MPSSR_mar2021_hg38/Pipeline_output/SeuratObject.rds") 

 

 

# Read in and add Metadata 

ctc_info <- read.table(file = "C:/Users/changmat/Documents/SeuratData/CTC_patient_info3_MPSSR_updated_patient_names.csv",row.names = NULL,header=T,sep=",") 

sample_meta <- ctc_info[5,] 

sample_meta <- sample_meta[2:97] 

sample_meta <- t(sample_meta) 

rownames(sample_meta) <- sort(colnames(x=SO)) 

SO <- AddMetaData(object= SO, metadata = sample_meta, col.name = "sample3") 

 

patient_status <- ctc_info[6,] 

patient_status <- patient_status[2:97] 

patient_status <- t(patient_status) 

rownames(patient_status) <- sort(colnames(x=SO)) 

SO <- AddMetaData(object= SO, metadata = patient_status, col.name = "patient_status") 

 

patient_id <- ctc_info[4,] 

patient_id <- patient_id[2:97] 

patient_id <- t(patient_id) 

rownames(patient_id) <- sort(colnames(x=SO)) 

SO <- AddMetaData(object=SO, metadata=patient_id, col.name="patient_id") 

 

names <- colnames(x=SO) 

SO <- AddMetaData(object=SO, metadata=names, col.name="names") 

 

 

#Subset the dataset based on >250 features expressed per cell and < 20% mitochondrial DNA per cell 

QC_subset <- subset(SO, subset = nFeature_RNA > 250 & percent.mt < 20) 

 

#Create UMAP plots by patient ID and patient status 

Idents(QC_subset) <- "patient_id" 

DimPlot(QC_subset)+ xlim(-7,7) + ylim(-3,3) 

Idents(QC_subset) <- "patient_status" 

DimPlot(QC_subset)+ xlim(-7,7) + ylim(-3,3) 

 

#Create QC metrics violin plot 

VlnPlot(QC_subset, features=c("nFeature_RNA", "percent.mt")) 

 

#Find differential expression based on patient status metadata subsets 

Idents(QC_subset) <- "patient_status" 

 

#Find differential expression of Patients vs PBMC clusters 

QC_patient_vs_PBMC_markers <- FindMarkers(QC_subset, group.by="patient_status",ident.1=c("Patient"), ident.2=c("PBMC"), only.pos = FALSE, min.pct = 0, log2fc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) 

 

#Find differential expression of Culture (A375 and Mini Bulk Prep) vs PBMC clusters 

QC_culture_vs_PBMC_markers <- FindMarkers(QC_subset, group.by="patient_status",ident.1=c("A375", "Mini_Bulk_Prep"), ident.2= "PBMC", only.pos = FALSE, min.pct = 0, log2fc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) 

 

#Find differential expression of Patient vs Culture clusters  

QC_combined_vs_PBMC_markers <- FindMarkers(QC_subset, group.by="patient_status",ident.2=c("Patient", "A375", "Mini_Bulk_Prep"), ident.1="PBMC", only.pos = FALSE, min.pct = 0, log2fc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) 

 

#Find differential expression of all samples vs PBMC clusters 

QC_pat_vs_culture_markers <- FindMarkers(QC_subset, group.by="patient_status",ident.1=c("Patient", "A375", "Mini_Bulk_Prep"), ident.2="PBMC", only.pos = FALSE, min.pct = 0, log2fc.threshold = 0, min.cells.feature = 0, min.cells.group = 0) 

 

#Take the top 10 most significant DE genes from each category 

QC_patient_vs_PBMC_markers_top10 <- QC_patient_vs_PBMC_markers %>% top_n(n = 10, wt = avg_log2FC) 

QC_culture_vs_PBMC_markers_top10 <- QC_culture_vs_PBMC %>% top_n(n = 10, wt = avg_log2FC) 

QC_combined_vs_PBMC_markers_top10 <- QC_combined_vs_PBMC_markers %>% top_n(n = 10, wt = avg_log2FC) 

QC_patient_vs_culture_markers_top10 <- QC_pat_vs_culture_markers %>% top_n(n = 10, wt = avg_log2FC) 

 

#Combine genes for a combined analysis heatmap 

comparison_genes <- c(rownames(QC_patient_vs_PBMC_markers_top10), rownames(QC_culture_vs_PBMC_markers_top10), rownames(QC_combined_vs_PBMC_markers_top10), rownames(QC_patient_vs_culture_markers_top10)) 

 

#Generate individual heatmaps for each comparison 

Idents(QC_subset) <- "patient_status" 

DoHeatmap(QC_subset, features=comparison_genes) 

DoHeatmap(QC_subset, features=rownames(QC_patient_vs_PBMC_markers_top10)) 

DoHeatmap(QC_subset, features=rownames(QC_culture_vs_PBMC_markers_top10)) 

DoHeatmap(QC_subset, features=rownames(QC_combined_vs_PBMC_markers_top10)) 

DoHeatmap(QC_subset, features=rownames(QC_patient_vs_culture_markers_top10)) 

 

#Generate feature plots for individual genes in the QC subset 

FeaturePlot(QC_subset, features="PTPRC")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="PRAME")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="CSPG4")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="MCAM")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="KIT")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="CD274")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="MLANA")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="MITF")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="CTNNB1")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="WNT5A")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="VIM")+ xlim(-7,7) + ylim(-3,3) 

FeaturePlot(QC_subset, features="VEGFA")+ xlim(-7,7) + ylim(-3,3) 

 

#Calculate statistics for nFeature and percent.mt for each group in patient_status metadata (Patient, A375, PBMC, Mini Bulk Prep) 

 

QC_A375_subset <- subset(QC_subset, idents = c("A375")) 

QC_PBMC_subset <- subset(QC_subset, idents = c("PBMC")) 

QC_Mini_subset <- subset(QC_subset, idents = c("Mini_Bulk_Prep")) 

QC_Patient_subset <- subset(QC_subset, idents = c("Patient")) 

 

sd(QC_Patient_subset@meta.data$nFeature_RNA) 

mean(QC_Patient_subset@meta.data$nFeature_RNA) 

 

sd(QC_A375_subset@meta.data$nFeature_RNA) 

mean(QC_A375_subset@meta.data$nFeature_RNA) 

 

sd(QC_PBMC_subset@meta.data$nFeature_RNA) 

mean(QC_PBMC_subset@meta.data$nFeature_RNA) 

 

sd(QC_Mini_subset@meta.data$nFeature_RNA) 

mean(QC_Mini_subset@meta.data$nFeature_RNA) 

 

sd(QC_Patient_subset@meta.data$percent.mt) 

mean(QC_Patient_subset@meta.data$percent.mt) 

 

sd(QC_A375_subset@meta.data$percent.mt) 

mean(QC_A375_subset@meta.data$percent.mt) 

 

sd(QC_PBMC_subset@meta.data$percent.mt) 

mean(QC_PBMC_subset@meta.data$percent.mt) 

 

sd(QC_Mini_subset@meta.data$percent.mt) 

mean(QC_Mini_subset@meta.data$percent.mt) 

 

 

#analysis for figure 7a 

 

#Subset to patient cells only 

Idents(QC_subset) <- "patient_status" 

QC_patient_subset <- subset(QC_subset, idents = c("Patient")) 

 

#Process patient cells separately 

QC_patient_subset <- NormalizeData(QC_patient_subset) 

QC_patient_subset <- FindVariableFeatures(QC_patient_subset, selection.method = "vst", nfeatures = 2000) 

QC_patient_subset <- ScaleData(QC_patient_subset) 

QC_patient_subset <- RunPCA(QC_patient_subset, npcs = 30, verbose = TRUE) 

QC_patient_subset <- RunUMAP(QC_patient_subset, reduction = "pca", dims = 1:9)  

QC_patient_subset <- FindNeighbors(QC_patient_subset) 

QC_patient_subset <- FindClusters(QC_patient_subset, resolution = 1.2) 

 

#Plot Patient Cells  

Idents(QC_patient_subset) <- QC_patient_subset@meta.data$patient_id 

DimPlot(QC_patient_subset, reduction = "pca") 

DimPlot(QC_patient_subset, reduction = "umap") 

 

#sPLS-DA analysis 

library(mixOmics) 

 

X <- t(QC_patient_subset@assays$RNA@data) 

X <- as.matrix(X) 

Y <- QC_patient_subset@meta.data$patient_id 

Y <- as.factor(Y) 

MyResult.splsda <- splsda(X,Y) 

 

#sPLS-DA Plot 

plotIndiv(MyResult.splsda, ind.names = F, legend=TRUE, 

          ellipse = T, star = T, centroid = T, title = 'sPLS-DA', 

          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2', alpha = 1, style = "lattice", pch = 15:22)  

 

loadings <- MyResult.splsda$loadings 

pat_loadings <- as.data.frame(loadings$Y) 

pat_indicator_mat <- as.data.frame(MyResult.splsda$ind.mat) 

pat_variate_mat <- as.data.frame(MyResult.splsda$variates) 

 

#Generate ANOSIM from sPLS-DA plot values 

anosim_vals <- anosim(pat_variate_mat, Y, distance = "bray") 

 

#Generate heatmap of sPLS-DA distances from centroids 

d <- dist(pat_variate_mat) 

d <- as.matrix(d) 

colnames(d) <- Y 

rownames(d) <- Y 

d <- d[, sort(colnames(d))] 

heatmap(d, keep.dendro = FALSE, verbose = F, scale = "row", col=heat.colors(5))# + 

plot.new() 

legend(x="bottom", legend=c(100, 90, 80,70, 60, 50, 40, 30, 20, 10, 0), fill=heat.colors(11), cex = 0.6) 

heatmap(d,Rowv=NA,Colv=NA,col=heat.colors(15)) 

 

#Sort sPLS-DA distances by Patient ID 

Pat_1_vals <- "" 

Pat_2_vals <- "" 

Pat_3_vals <- "" 

Pat_4_vals <- "" 

Pat_5_vals <- "" 

Pat_6_vals <- "" 

Pat_7_vals <- "" 

 

for(i in 1:length(colnames(d))){ 

  if(colnames(d)[i] == levels(factor(colnames(d)))[1]){ 

    Pat_1_vals <- c(Pat_1_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[2]){ 

    Pat_2_vals <- c(Pat_2_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[3]){ 

    Pat_3_vals <- c(Pat_3_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[4]){ 

    Pat_4_vals <- c(Pat_4_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[5]){ 

    Pat_5_vals <- c(Pat_5_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[6]){ 

    Pat_6_vals <- c(Pat_6_vals, d[,i]) 

  } 

  else if(colnames(d)[i] == levels(factor(colnames(d)))[7]){ 

    Pat_7_vals <- c(Pat_7_vals, d[,i]) 

  } 

} 

 

Pat_1_vals <- Pat_1_vals[-1] 

Pat_2_vals <- Pat_2_vals[-1] 

Pat_3_vals <- Pat_3_vals[-1] 

Pat_4_vals <- Pat_4_vals[-1] 

Pat_5_vals <- Pat_5_vals[-1] 

Pat_6_vals <- Pat_6_vals[-1] 

Pat_7_vals <- Pat_7_vals[-1] 

 

 

#Boxplot of distances from centroid on sPLS-DA plot by Patient ID 

boxplot(as.numeric(Pat_1_vals), as.numeric(Pat_2_vals), as.numeric(Pat_3_vals), as.numeric(Pat_4_vals), as.numeric(Pat_5_vals), as.numeric(Pat_6_vals), as.numeric(Pat_7_vals), main = "Euclidian Distance from Centroid by Patient", names = c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5", "Patient 6", "Patient 7")) 

 

#Generate data frame to hold values by Patient ID 

vals_by_pat <- c(Pat_1_vals, Pat_2_vals, Pat_3_vals, Pat_4_vals, Pat_5_vals, Pat_6_vals, Pat_7_vals) 

 

vals_by_pat <- as.data.frame(vals_by_pat) 

Patient <- "" 

Patient[1:2025] <- "Patient 1" 

Patient[2026:3236] <- "Patient 2" 

Patient[3236:3311] <- "Patient 3" 

Patient[3311:3761] <- "Patient 4" 

Patient[3761:3911] <- "Patient 5" 

Patient[3911:5561] <- "Patient 6" 

Patient[5561:5625] <- "Patient 7" 

vals_by_pat$Patient <- Patient 

 

#Tukey Test on distance values 

res.aov <- aov(Values ~ Patient, data = vals_by_pat) 

summary(res.aov) 

TukeyHSD(res.aov) 

 
