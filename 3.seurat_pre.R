library(Seurat)
library(harmony)
set.seed(123)

sce_list <- list()

sce_list[[sample_id]] <- readRDS("finddoublet.rds")

sce <- merge(sce_list[[1]],sce_list[-1],add.cell.ids =names(sce_list))

sce_filter <- seurat_preprocess(sce_filter,sct = FALSE,var.to.regress = TRUE,vars = c("percent.mt","nCount_RNA","nFeature_RNA"),scale.all.gene = TRUE)

sce_filter <- RunHarmony(sce_filter,group.by.vars = "orig.ident",reduction.save = "harmony")

pc_num_select <- pc_select(sce_filter,ndim = 50,reduction = "harmony")

# pc_num_select <- 20
sce_filter <- find_neighbors_cluster_umap(sce_filter,sct = F,algorithm = 1,reduction = "harmony",cluster.resns=c(seq(0,1,0.1)),ndim = pc_num_select,umap.reduction.name = "harmony.umap")

sce_filter <- JoinLayers(sce_filter)
sce_filter_all_marker <- FindAllMarkers(sce_filter,only.pos = T,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25) %>% mutate(difference = .$pct.1 - .$pct.2)
sce_filter_all_marker_deg <- sce_filter_all_marker %>% group_by(cluster) %>% top_n(n = 250,wt = avg_log2FC) %>% arrange(cluster,desc(avg_log2FC))

write_csv(sce_filter_all_marker_deg,file = "deg.csv")