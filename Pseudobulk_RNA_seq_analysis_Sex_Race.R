# Load libraries  R version 4.1.3
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
#install.packages("cowplot")
#install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
library(Matrix.utils)
library(edgeR)
library(dplyr)
#install.packages("magrittr")
library(magrittr)
#install.packages("Matrix")
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
#BiocManager::install("apeglm")
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Bring in Seurat object
seurat <- readRDS("Donor_NICM_Final.rds")
table(seurat$orig.ident)
#################remoce the underscore of orig.ident or It will cause error for smaple_id extract
seurat$orig.ident<-gsub("Donor_","Donor",seurat$orig.ident)
seurat$orig.ident<-gsub("DCM_","DCM",seurat$orig.ident)
table(seurat$orig.ident)
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 
metadata <- seurat@meta.data
Idents(seurat)<-seurat$celltype
seurat@active.ident
# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)
metadata$sample_id<-factor(seurat$orig.ident)
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Explore the raw counts for the dataset
## Check the assays present
assays(sce)
## Explore the raw counts for the dataset
dim(counts(sce))
counts(sce)[1:6, 1:6]
## Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))
#First, we need to determine the number of clusters and the cluster names present in our dataset.
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns
# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei
## Perform QC if not already performed
#dim(sce)

# Calculate quality control (QC) metrics
sce2 <- calculateQCMetrics(sce)

# Get cells w/ few/many detected genes
#sce2$is_outlier <- isOutlier(
 # metric = sce2$total_features_by_counts,
 # nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
#sce2 <- sce2[, !sce2$is_outlier]
#dim(sce2)

## Remove lowly expressed genes which have less than 10 cells with any counts
#sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
#dim(sce2)

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

splitf
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)
# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)
# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()
# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()
# Create a data frame with the sample IDs, cluster IDs and condition(group or Sex or Race)

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group) 

metadata       
# Generate vector of cluster IDs
class(metadata)
stringsAsFactors=T
factor(metadata$cluster_id)
#lapply(metadata, levels)
#levels(metadata[,1])
#clusters <- table(metadata[,1])
clusters<-levels(factor(metadata[,1]))
#clusters<-c("adipocytes","Cardiomyocytes",
#"Endo","Fibro","Gila","lymphocytic vessels","Mast Cell","Myeloid","Pericyte&SMC","NK&T")
clusters
save(pb,metadata,counts,ei,sce,samples_list,groups,file = "Donor_NICMCM_pseudotime.Rdata")
load("Pesudobulk-analysis.Rdata")

for (i in 1:length(clusters)) {
  
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[i]), ]
  head(cluster_metadata)
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  head(cluster_metadata)
  
  # Subset the counts to only the cardiomyocytes
  
  
  counts <- pb[[clusters[i]]]
  
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  all(rownames(cluster_metadata) == colnames(cluster_counts))        
  
  # keep all celltype for DEseq analysis
  
  #creat DEseq object 
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group)
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  
  # Plot PCA
  
  DESeq2::plotPCA(rld, intgroup = "group")
  ggsave(paste0(clusters[i],"based_pcaplot.png"),width = 5,height=5,path = "output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/")
  dev.off()
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  # Plot heatmap
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop=F])
  ggsave(paste0(clusters[i],"based_heatmap.png"),width = 5,height=5,path = "output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/")
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  # Plot dispersion estimates
  #plotDispEsts(dds)
  #ggsave(paste0(clusters[i],"based_dispersion estimates.png"),width = 5,height=5,path = "output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/")
  
  # Output results of Wald test for contrast for Dox vs ctrl
  levels(as.factor(cluster_metadata$group))[2]
  levels(as.factor(cluster_metadata$group))[1]
  
  contrast <- c("group", levels(as.factor(cluster_metadata$group))[2], levels(as.factor(cluster_metadata$group))[1])
  contrast
  # resultsNames(dds)
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05)
  resultsNames(res) 
  #BiocManager::install("ashr")
  library(ashr)
  res <- lfcShrink(dds, 
                   contrast =  contrast,res = res,type="ashr")
  #First let’s generate the results table for all of our results:
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # Check results output
  res_tbl
  
  # Write all results to file
  write.csv(res_tbl,
            paste0("output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/", clusters[i], "_", levels(as.factor(cluster_metadata$group))[2], "_vs_", levels(as.factor(cluster_metadata$group))[1], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Set thresholds
  padj_cutoff <- 0.05
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  # Check significant genes output
  sig_res
  
  # Write significant results to file
  write.csv(sig_res,
            paste0("output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/", clusters[i], "_", levels(as.factor(cluster_metadata$group))[2], "_vs_", levels(as.factor(cluster_metadata$group))[1], "_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  ## ggplot of top genes
  normalized_counts <- counts(dds, 
                              normalized = TRUE)
  
  ## Order results by padj values
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n=20)
  
  
  top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_genes)
  
  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
  
  gathered_top20_sig <- inner_join(ei[, c("sample_id", "group" )], gathered_top20_sig, by = c("sample_id" = "samplename"))
  
  ## plot using ggplot2
  p<-ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene, 
                   y = normalized_counts, 
                   color = group), 
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  ggsave(paste0(clusters[i],"based_Top20.pdf"),width = 10,height=5,path = "output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/")
  
  # Extract normalized counts for only the significant genes
  sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% sig_res$gene)
  
  # Set a color palette
  heat_colors <- brewer.pal(6, "YlOrRd")
  
  # Run pheatmap using the metadata data frame for the annotation
  # pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
  #          color = heat_colors, 
  #          cluster_rows = T, 
  #          show_rownames = F,
  #          annotation = cluster_metadata[, c("group", "cluster_id")], 
  #          border_color = NA, 
  #           fontsize = 10, 
  #           scale = "row", 
  #           fontsize_row = 10, 
  #           height = 20)        
  
  ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
  res_table_thres <- res_tbl %>% 
    mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  
  ## Volcano plot
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    ggtitle("Volcano plot of  relative to group") +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))      
  ggsave(paste0(clusters[i],"volconoplot.pdf"),width = 10,height=5,path = "output/RNA_longnormilization/pseudobulk rna-seq/Donor_NICM/")
  save(dds,pb,ei,metadata,clusters,kids,sids,splitf,file = paste0(clusters[i],"_group_Pesudobulk-analysis.Rdata"))
  
}





