args <- commandArgs()

help <- function(){
  cat("accessibility_profiles.R :
Make profiles plots around promoters. Output will be save in the parsed directory in the covPath. The promter and body regions (gene body w#ithout promoter) will be saved in the outDir.\n")
  cat("Usage: \n")
  cat("--covThresh   : Minimum library size                                        [required]\n")
  cat("--featThresh  : Minimum number of expressed features                        [required]\n")
  cat("--top50       : Maximum fraction of reads accounting for top 50 features    [required]\n")
  cat("--featureLQ   : Feature lower quantile                                      [required]\n")
  cat("--featureUQ   : Feature upper quantile                                      [required]\n")
  cat("--CountUQ     : Count upper quantile                                        [required]\n")
  cat("--MTUQ        : percent MT upper quantile                                   [required]\n")
  cat("--integrate   : TRUE or FALSE determining whether to integrate              [required]\n")
  cat("--res         : resolution for clustering                                   [required]\n")
  cat("\n")
  q()
}

io   <- list()
opts <- list()

# Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
  help()
} else {
  coverage_threshold      <- sub( '--covThresh=', '', args[grep('--covThresh=', args)] )
  features_threshold      <- sub( '--featThresh=', '', args[grep('--featThresh=', args)] )
  top50_threshold         <- sub( '--top50=', '', args[grep('--top50=', args)] )
  Feature_lowerQuantile   <- sub( '--featureLQ=', '', args[grep('--featureLQ=', args)])
  Feature_upperQuantile   <- sub( '--featureUQ=', '', args[grep('--featureUQ=', args)])
  Count_upperQuantile     <- sub( '--CountUQ=', '', args[grep('--CountUQ=', args)])
  percentMT_upperQuantile <- sub( '--MTUQ=', '',args[grep('--MTUQ=',args)])
  integrateTF             <- sub( '--integrate=', '',args[grep('--integrate=',args)])
  res                     <- sub( '--res=', '', args[grep('--res=',args)])
  
}


opts <- list()
opts$coverage_threshold <- as.numeric(coverage_threshold)

opts$features_threshold <- as.numeric(features_threshold)

opts$top50_threshold <- as.numeric(top50_threshold)

opts$Feature_lowerQuantile <- as.numeric(Feature_lowerQuantile)

opts$Feature_upperQuantile <- as.numeric(Feature_upperQuantile)

opts$Count_upperQuantile <- as.numeric(Count_upperQuantile)

opts$percentMT_upperQuantile <- as.numeric(percentMT_upperQuantile)

opts$integrate <- integrateTF

opts$res <- as.numeric(res)

## I/O ##
io <- list()
io$in.gene_metadata   <- "data/gene_metadata.tsv"
io$in.sample_metadata <- "data/counts/sample_metadata.tsv"
io$in.raw_counts      <- "data/counts/raw_counts_.filt.tsv"
io$out.file           <- "data/seurat/SeuratObject.rds"
io$out.file.CCred     <- "data/seurat/SeuratObject_CCred.rds"
io$dataDir            <- "data/seurat"
io$plotDir            <- "plots/seurat"

if( !( file.exists(io$dataDir)) )  {
  dir.create( io$dataDir, FALSE, TRUE )  
}

if(!( file.exists( io$plotDir ) ) ) {
  dir.create( io$plotDir, FALSE, TRUE )  
}

library(Seurat)
library(data.table)
library(purrr)
library(scater)
library(ggplot2)
library(cowplot)
library(png)
library(dplyr)
library(stringr)


fread_df <- partial(fread, data.table = FALSE)

## read in counts data and associated metadata ##
counts <- fread_df(io$in.raw_counts) %>% 
  tibble::column_to_rownames("ens_id")
counts$Genes <- NULL

feature_metadata           <- fread_df(io$in.gene_metadata)
rownames(feature_metadata) <- feature_metadata$ens_id
stopifnot(rownames(counts) == feature_metadata$ens_id)

counts           <- as.matrix(counts)
rownames(counts) <- feature_metadata$gene

## Create Seurat Object ##
#rownames(counts) <- feature_metadata[rownames(counts),3]
SO <- CreateSeuratObject(counts = counts, min.cells = 0,
                         min.features = 0)

## save unique gene ids  
feature_metadata$gene_unique <- rownames(SO)

ctr <- 0
SO@meta.data$origin <- str_extract(rownames(SO@meta.data), "[^_]+")

SO[["percent.mt"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(paste0(io$plotDir, "/QC_vlnplot", as.character(ctr), ".pdf"), plot1)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot2 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
save_plot(paste0(io$plotDir, "/QC_scatter_percentMT", as.character(ctr), ".pdf"), plot2)

plot3 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
save_plot(paste0(io$plotDir, "/QC_scatter_nfeature", as.character(ctr), ".pdf"), plot3)


SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features
SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SO), 10)

plot4 <- VariableFeaturePlot(SO)
#plot5 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(paste0(io$plotDir, "/VariableFeaturePlot", as.character(ctr), ".pdf"), plot4)

#Scaling the data
all.genes <- rownames(SO)
SO <- ScaleData(SO, features = all.genes)

#Perform linear dimensional reduction

SO <- RunPCA(SO, features = VariableFeatures(object = SO))

plot6 <- DimHeatmap(SO, dims = 1:6, cells = 500, balanced = TRUE)
save_plot(paste0(io$plotDir, "/Heatmap1", as.character(ctr), ".pdf"), plot6)


SO <- JackStraw(SO, num.replicate = 100)
SO <- ScoreJackStraw(SO, dims = 1:20)
JackStrawPlot(SO, dims = 1:20)

plot7 <- ElbowPlot(SO)
save_plot(paste0(io$plotDir, "/Elbowplot", as.character(ctr), ".pdf"), plot7)

#Cluster the cells
SO <- FindNeighbors(SO, dims = 1:9)

SO <- FindClusters(SO, resolution = 0.2)

#Run non-linear dimensional reduction
#UMAP

SO <- RunUMAP(SO, dims = 1:10)


plot8 <- DimPlot(SO, reduction = "umap", pt.size = 4)
save_plot(paste0(io$plotDir, "/UMAP", as.character(ctr), ".pdf"), plot8)

saveRDS(SO, io$out.file)

SO.markers <- FindAllMarkers(SO, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
SO.top10.markers <- SO.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


plot9 <- DoHeatmap(SO, features = SO.top10.markers$gene)
save_plot(paste0(io$plotDir, "/top10_heatmap", as.character(ctr), ".pdf"), plot9)

plot10 <- VlnPlot(SO, 
        features = SO.top10.markers$gene, 
        stack = TRUE, 
        same.y.lims = FALSE, 
        pt.size = 0, 
        flip = TRUE)
save_plot(paste0(io$plotDir, "/top10_vlnplot", as.character(ctr), ".pdf"), plot10)



