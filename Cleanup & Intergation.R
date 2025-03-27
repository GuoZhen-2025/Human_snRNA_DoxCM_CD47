#### 0. Set pathway ####
setwd("//Project")

#### 1. Jupyter run scrublet to calculate doublet score sample by sample ####
# Perform in python "Individual Doublet Removal.ipynb"

#### 2. library packages ####
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)
library(harmony)
library(ggstatsplot)
library(tidyverse)
library(Nebulosa)

#### 3. Unzip 10x genomics data for doublet prediction in python ####
gunzip("./XYZ/matrix.mtx.gz", remove = F)
gunzip("./XYZ/features.tsv.gz", remove = F)
gunzip("./XYZ/barcodes.tsv.gz", remove = F)

#### 4. Create seuratobj for individual ####
# Sample_xxxx
gc()
dir.create("./Sample_xxxx")
sample_dir <- "./Sample_xxxx/"
sample.data <- Read10X(data.dir =sample_dir)
sample <- CreateSeuratObject(counts = sample.data)
# Add group info and calculate mitochondria
sample$orig.ident <- "xxxx"
sample$group <- "DoxCM" # Donor or NICM
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
# Read doublet prediction file
pred <- read.csv("./Sample_xxxx/unzip/doublet_prediction.txt", header = T, sep = ",")
pred[pred$prediction == "False",]$prediction <- "Singlet"
pred[pred$prediction == "True",]$prediction <- "Doublet"
table(pred$prediction)
# Remove doubletes
sample@meta.data$scrublet <- pred$prediction
final_sample <- subset(sample, scrublet == "Singlet")
# Save final RDS
saveRDS(final_sample, "./Sample_xxxx/Final_Sample_xxxx_removal_scrublet.rds")

#### 5. Merge all the seuratobj ####
# load files 
rm = ls()
sample01 <- readRDS("./Sample_xxxx/Final_Sample_xxxx_removal_scrublet.rds")
...
sample25 <- readRDS("./Sample_xxxx/Final_Sample_xxxx_removal_scrublet.rds")
# marge
merged_seuratobj <- merge(sample01, y = c(sample02, sample03, ..., sample25))
saveRDS(merged_seuratobj, file = "./DoxCM_3groups_rawseuratobj.rds") ### rawdata named as "DoxCM_3groups_rawseuratobj.rds" --- rawdata 1

#### 6. QC ####
seuratobj <- merged_seuratobj
VlnPlot(seuratobj, features = c("nFeature_RNA"), pt.size = 0.0001, group.by = "orig.ident")
VlnPlot(seuratobj, features = c("nCount_RNA"), pt.size = 0.0001, group.by = "orig.ident")
VlnPlot(seuratobj, features = c("percent.mt"), pt.size = 0.0001, group.by = "orig.ident")
seuratobj$log10GenesPerUMI <- log10(seuratobj$nFeature_RNA)/log10(seuratobj$nCount_RNA)
VlnPlot(seuratobj, features = c("log10GenesPerUMI"), pt.size = 0.0001, group.by = "orig.ident")

Fseuratobj <- subset(seuratobj, subset = nCount_RNA <= 40000 &
                       nFeature_RNA >= 250 &  nFeature_RNA < 10000 &
                       log10GenesPerUMI > 0.80 &
                       percent.mt < 8)
VlnPlot(Fseuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), pt.size = 0, ncol = 2, group.by = "orig.ident")

#### 7. Normalization ####
Fseuratobj <- NormalizeData(Fseuratobj)
Fseuratobj <- FindVariableFeatures(Fseuratobj)
Fseuratobj <- ScaleData(Fseuratobj)
Fseuratobj <- RunPCA(Fseuratobj)
ElbowPlot(Fseuratobj, ndims = 50)

Fseuratobj <- FindNeighbors(Fseuratobj, dims = 1:50)
Fseuratobj <- FindClusters(Fseuratobj)
Fseuratobj <- RunUMAP(Fseuratobj, dims = 1:50)

#### 8. Intergation with Harmonoy ####
### Harmony
Fseuratobj <- readRDS("./DoxCM_3groups_SCT.rds")
DefaultAssay(Fseuratobj) <- "RNA"
Fseuratobj <- SCTransform(Fseuratobj, vars.to.regress = c("percent.mt","group"))

Fseuratobj <- RunHarmony(Fseuratobj, group.by.vars = c("group","orig.ident"), assay.use = "SCT")
Fseuratobj <- RunUMAP(Fseuratobj, reduction = "harmony", dims = 1:50)
Fseuratobj <- FindNeighbors(Fseuratobj, reduction = "harmony", dims = 1:50) %>% FindClusters()

### Cleanup: using CellSelector and subset functions in 11 separated major clusters to exclude junk and high mitochondrial cells
### Merge and rerun step 8 get the rawdata named as "DoxCM_cleanup_combined_SCT_harmony_by_clusters.rds" --- rawdata 2

