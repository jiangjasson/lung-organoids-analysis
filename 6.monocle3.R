library(monocle3)

cds <- sce_filter
data <- GetAssayData(cds,assay = "RNA",layer = "counts")
cell_metadata <- cds@meta.data#[,c(1,2,3,5,6,7,8,9,10,11,13,33,34,52)]
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

cds <- preprocess_cds(cds,num_dim = 100) 
plot_pc_variance_explained(cds) #elbowplot
cds <- reduce_dimension(cds,preprocess_method = "PCA",reduction_method = "UMAP")
cds <- cluster_cells(cds)

#seurat umap
##1.assign partitions
recreate.partition <- rep(1,length(cds@colData@rownames)) 
names(recreate.partition) <- cds@colData@rownames 
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- sce_filter$celltype
cds@clusters$UMAP$clusters <- list_cluster

##3. seurat的UMAP坐标信息分配给cds----
cds@int_colData@listData$reducedDims$UMAP <- sce_filter@reductions$harmony.umap@cell.embeddings

p <- plot_cells(cds,color_cells_by = "partition",
                label_branch_points = F,
                label_leaves = F,
                cell_size = 1,
                label_groups_by_cluster = F,group_label_size = 3,graph_label_size = 3
)

cds <- learn_graph(cds)
p1 <- plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = F,
                 label_leaves = F,
                 label_branch_points = T,
                 group_label_size = 4,graph_label_size = 3,
                 cell_size = 0.5)

cds <- order_cells(cds = cds)
p2 <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups = F,
                label_leaves = T,
                label_branch_points = T,
                graph_label_size = 3,
                group_label_size = 3,cell_size = 1.5)

pseudotime <- pseudotime(cds)
# sce_filter$pseudotime <- 0
sce_filter$pseudotime <- pseudotime[colnames(sce_filter)]