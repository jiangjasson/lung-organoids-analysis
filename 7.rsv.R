library(Seurat)
library(harmony)
set.seed(123)

rsv_gene <- c("NS1","NS2","N","P","M","SH","G","F","M2","L")

sce_list <- list()

sce_list[[sample_id]] <- readRDS("finddoublet.rds")

for (i in 1:length(sce_list)){
  sce_list[[i]]$percent.rsv <- PercentageFeatureSet(sce_list[[i]],features = rsv_gene)
}

sce <- merge(sce_list[[1]],sce_list[-1],add.cell.ids =names(sce_list))

sce_filter <- seurat_preprocess(sce_filter,sct = FALSE,var.to.regress = TRUE,vars = c("percent.mt","nCount_RNA","nFeature_RNA"),scale.all.gene = TRUE)

sce_filter <- find_neighbors_cluster_umap(sce_filter,sct = F,algorithm = 1,reduction = "pca",cluster.resns=c(seq(0,1,0.1)),ndim = pc_num_select,umap.reduction.name = "umap")
