library(Seurat)
library(DoubletFinder)
set.seed(123)

#load
wd <- getwd()
sample_dir <- list.files(path = paste0(wd,"/1.Matrix"),full.names = TRUE)
sample_list <- list()
for (i in sample_dir){
  counts <- Read10X(i)
  sample_id <- basename(i)
  seurat_obj <- CreateSeuratObject(counts = counts,project = sample_id)
  sample_list[[sample_id]] <- seurat_obj
}

sce_list <- sample_list
for (i in 1:length(sce_list)){
  sce_list[[i]]$percent.mt <- PercentageFeatureSet(sce_list[[i]],pattern = "^MT-")
  sce_list[[i]]$percent.ribo <-PercentageFeatureSet(sce_list[[i]],pattern = "^RP[SL]")
  sce_list[[i]]$percent.hb <- PercentageFeatureSet(sce_list[[i]],pattern = "^HB[^(P)]")
}

sce_filter <- lapply(sce_list,FUN = function(seu){
  seu <- seurat_filter_cell(seu,nfeature_lower = 300,nfeature_quantile = TRUE,percentmt = 10,keep_quantile = 0.95,merge = TRUE)
})

sce_list <- lapply(sce_list,FUN = function(seu){
  seu <- seurat_filter_gene(seu,merge = FALSE,ngene = 3,multi_sample_counts_in_layers = FALSE)
})
  
sce_filter <- lapply(sce_list,FUN = function(seu){
  seu <- seurat_preprocess(seu,sct = FALSE,var.to.regress = TRUE,vars = c("percent.mt","nFeature_RNA","nCount_RNA"),scale.all.gene = FALSE)
})

sce_list <- lapply(sce_list,FUN = function(seu){
  pc_num_select <- pc_select(seu,reduction = "pca")
  #pc_num_select <- 15
  seu <- find_neighbors_cluster_umap(seu,reduction = "pca",sct = FALSE,cluster.resns = 0.5,ndim = pc_num_select,algorithm = 1,umap.reduction.name = "umap")
})

sce_list <- lapply(sce_list,FUN = function(x){
  #x@meta.data$seurat_clusters <- x@meta.data$SCT_snn_res.0.5
  x <- run_doublefinder(x,sct = FALSE)
})

for (id in seq_along(sce_list)){
  sce_list[[id]]@meta.data$Doublet <- "NA"
  sce_list[[id]]@meta.data[,"Doublet"] <- sce_list[[id]]@meta.data[,10]
  sce_list[[id]]@meta.data$Doublet[which(sce_list[[id]]@meta.data$Doublet == "Doublet" & sce_list[[id]]@meta.data[,11] == "Singlet")] <- "Doublet_low"
  sce_list[[id]]@meta.data$Doublet[which(sce_list[[id]]@meta.data$Doublet == "Doublet" )] <- "Doublet_high"
}

#提取矩阵，保存为单个文件
for (i in names(sce_list)){
  seu <- GetAssayData(sce_list[[i]],assay = "RNA",layer = "counts")
  seu <- CreateSeuratObject(seu,meta.data = sce_list[[i]]@meta.data["Doublet"],project = i)
  saveRDS(seu,file = paste0("../1.Matrix/",i,"finddoublet.rds"))
}
  