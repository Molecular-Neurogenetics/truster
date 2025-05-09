#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(patchwork)
set.seed(10)

option_list = list(
  make_option(c("-i", "--inpath"), type="character", default=NULL,
              help="RData paths", metavar="character"),
  make_option(c("-s", "--ids"), type="character", default=NULL,
              help="Sample ids in order of appearance in -i", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path", metavar="character"),
  make_option(c("-e", "--experimentName"), type="character", default=NULL,
              help="Experiment name", metavar="character"),
  make_option(c("-r", "--resolution"), type="character", default=0.5,
              help="Resolution. Default 0.5", metavar="character"),
  make_option(c("-n", "--normalizationMethod"), type="character", default="LogNormalize",
              help = "Seurat normalization method (LogNormalize | CLR | RC)", metavar = "character"),
  make_option(c("-I", "--integrateSamples"), type="character", default = "FALSE", 
              help = "Integrate samples (Integration of multiple datasets, correct for batch effects or technology differences) (TRUE | FALSE). Default: FALSE", metavar = "character"),
  make_option(c("-S", "--maxSize"), type="character", default=500,
              help = "Maximum size for global variables in MiB (future.globals.maxSize). This will increase your RAM usage.", metavar = "numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) ){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

options(future.globals.maxSize = as.numeric(as.character(opt$maxSize)) * 1024^2)
paths <- trimws(unlist(str_split(opt$inpath, ',')))
ids <- trimws(unlist(str_split(opt$ids, ',')))
outpath <- ifelse(endsWith(opt$outpath, "/"), opt$outpath, paste(opt$outpath, '/', sep=''))
experiment_name <- trimws(opt$experimentName)
normalization_method <- as.character(opt$normalizationMethod)
integrate <- as.logical(opt$integrateSamples)
res <- as.numeric(as.character(opt$resolution))

print(c("Input paths: ", paths))
print(c("Input ids: ", ids))
print(c("Output path: ", outpath))
print(c("Experiment name: ", experiment_name))
print(c("Resolution: ", res))
print(c("Normalization method: ", normalization_method))
print(c("Integrate samples: ", integrate))

samples <- list()
for(i in 1:length(paths)){
  path <- paths[i]
  name <- ids[i]
  sample <- readRDS(path)
  sample <- AddMetaData(sample, Cells(sample), col.name = "original_cellIds")
  samples[[name]] <- sample
}

if(integrate == TRUE){
  experiment <- FindIntegrationAnchors(object.list = samples, dims = 1:30)
  experiment <- IntegrateData(anchorset = experiment, dims = 1:30)
  assay = "integrated"
}else{
  experiment <- merge(samples[[ids[1]]], y=samples[ids[-1]])
  assay = "RNA"

  experiment <- NormalizeData(experiment, normalization.method = normalization_method, assay = assay, margin = 2) 
  experiment <- FindVariableFeatures(experiment, selection.method = "vst", nfeatures = 2000, assay = assay)
}
rm(samples)
gc()
all.genes <- rownames(experiment)
experiment <- ScaleData(experiment, features = all.genes, assay = assay)
# experiment[["percent.mt"]] <- PercentageFeatureSet(experiment, pattern = "^MT-")

experiment <- RunPCA(experiment, features = VariableFeatures(object = experiment), assay = assay)
experiment <- FindNeighbors(experiment, dims = 1:10, assay = assay)
experiment <- FindClusters(experiment, resolution = res)
experiment <- RunUMAP(experiment, dims = 1:10, assay = assay)

colours <- colorRampPalette(brewer.pal(8, "Accent"))(length(unique(experiment$seurat_clusters)))
names(colours) <- as.character(unique(experiment$seurat_clusters))
experiment_colours <- merge(data.frame(seurat_clusters=unique(experiment$seurat_clusters)), data.frame(cluster_colours = colours, seurat_clusters = names(colours)), by='seurat_clusters')
seurat_clusters_df <- as.data.frame(experiment[["seurat_clusters"]])
seurat_clusters_df$Row.names <- rownames(seurat_clusters_df)
experiment_colours <- merge(seurat_clusters_df, experiment_colours, all.x = T)
rownames(experiment_colours) <- experiment_colours$Row.names
experiment_colours <- structure(as.character(experiment_colours$cluster_colours), names = as.character(experiment_colours$Row.names))
experiment <- AddMetaData(experiment, experiment_colours, col.name = "cluster_colours")

for(i in 1:length(ids)){
  id <- ids[i]
  sample <- subset(experiment, subset = orig.ident == id)
  embedding <- Embeddings(sample, reduction = "umap")
  
  embedding_origcellIds <- merge(embedding, experiment[['original_cellIds']][rownames(embedding),,drop=F], by='row.names')
  rownames(embedding_origcellIds) <- embedding_origcellIds$original_cellIds
  embedding_origcellIds <- embedding_origcellIds[,c("umap_1", "umap_2")]
  
  cluster_colours <- experiment[[c("original_cellIds", "seurat_clusters", "cluster_colours")]][rownames(embedding),,drop=F]
  rownames(cluster_colours) <- cluster_colours$original_cellIds
  cluster_colours <- cluster_colours[,-1]
  
  write.csv(embedding_origcellIds, file = paste(outpath, '/', id, "_cell_embeddings.csv", sep=''))
  write.csv(cluster_colours, file = paste(outpath, '/', id, "_clusters.csv", sep=''))
}

cluster_colours <- experiment[[c("original_cellIds", "seurat_clusters", "cluster_colours")]]
write.csv(cluster_colours, file = paste(outpath, '/', experiment_name, "_clusters.csv", sep=''))

embedding <- Embeddings(experiment, reduction = "umap")
write.csv(embedding, file = paste(outpath, '/', experiment_name, "_embeddings.csv", sep=''))

gene_expression <- FetchData(experiment, vars = rownames(experiment))
write.csv(gene_expression, file = paste(outpath, '/', experiment_name, "_gene_counts.csv", sep=''))

experiment[["original_cellIds"]] <- NULL
experiment[["cluster_colours"]] <- NULL

df <- as.data.frame(experiment$seurat_clusters)
colnames(df) <- 'clusters'


for(i in 1:length(ids)){
  sampleid <- ids[i]
  sample <- subset(experiment, subset = orig.ident == sampleid)
  
  df <- as.data.frame(sample$seurat_clusters)
  colnames(df) <- 'clusters'
  
  for (k in 1:length(unique(df$clusters))){
    cluster <- unique(df$clusters)[k]
    df.cluster <- subset(df, df$clusters == as.character(cluster))
    # Assuming the naming of the cells is [sample]_[barcode]
    df.cluster$barcode <- sapply(str_split(rownames(df.cluster), paste(unique(sample$orig.ident), "_", sep="")), `[[`, 2)
    
    dir.create(outpath)
    file_name <- paste(outpath, '/', sampleid, '_merged.clusters_', cluster, '.tsv', sep = '')
    print(file_name)
    # write.table(sapply(str_split(df.cluster$barcode, "_"), `[[`, 1), file = file_name, row.names = F, col.names = F, quote = F)
    write.table(df.cluster$barcode, file = file_name, row.names = F, col.names = F, quote = F)
  }
  
}

saveRDS(experiment, file=paste(outpath, experiment_name, ".rds", sep=''))