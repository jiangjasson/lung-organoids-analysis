#function
#doublefinder
run_doublefinder <- function(seu,sct = FALSE){
  ## pK 值确认 ---------------------------------------------------------------------------------------
  sweep.res.list_obj <- paramSweep(seu, PCs = 1:20, sct = sct)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  pK_value <- as.numeric(as.character(bcmvn_obj$pK[bcmvn_obj$BCmetric == max(bcmvn_obj$BCmetric)]))
  
  #同型双细胞评估----
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #找双细胞----
  pANN_value <- paste0("pANN_0.25_",pK_value,"_",nExp_poi)
  seu <- doubletFinder(seu, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  seu <- doubletFinder(seu, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = sct)
  return(seu)
}

seurat_filter_cell <- function(seu,nfeature_lower = 300,nfeature_uper = 10000,nfeature_quantile = FALSE,percentmt = 10,keep_quantile = 0.95,merge=FALSE){
  sample_name <- ifelse(merge, "",paste(unique(seu$orig.ident)))
  print(paste(sample_name,"Original cell count:", ncol(seu)))
  nfeature_uper <- ifelse(!nfeature_quantile,nfeature_uper,quantile(seu$nFeature_RNA, keep_quantile))
  seu <- subset(seu,subset = nFeature_RNA > nfeature_lower & 
                  nFeature_RNA < nfeature_uper &
                  percent.mt < percentmt & 
                  nCount_RNA < quantile(nCount_RNA, keep_quantile)
  )
  print(paste(sample_name,"Filtered cell count:", ncol(seu)))
  return(seu)
}

#过滤基因
seurat_filter_gene <- function(seu,ngene = 3,merge=FALSE,multi_sample_counts_in_layers = FALSE){
  sample_name <- ifelse(merge, "",paste(unique(seu$orig.ident)))
  if (multi_sample_counts_in_layers){
    seu <- JoinLayers(seu)
    print(paste("Original gene count:", nrow(seu)))
    total_counts <- GetAssayData(seu, layer="counts",assay = "RNA")
    # 计算每个基因在多少个细胞中表达
    gene_counts <- Matrix::rowSums(total_counts > 0)
    # 筛选出在至少n个细胞中表达的基因
    genes_to_keep <- names(gene_counts[gene_counts >= ngene])
    # 使用 subset 函数筛选出这些基因
    seu <- subset(seu, features = genes_to_keep)
    print(paste("Filtered gene count:", nrow(seu)))
    seu <- SplitObject(seu,split.by = "orig.ident")
    seu <- merge(x = seu[[1]],y = seu[-1],add.cell.ids =names(seu))
    return(seu)
  }else{
    print(paste(sample_name,"Original gene count:", nrow(seu)))
    keep_gene <- rowSums(seu@assays$RNA$counts > 0) > ngene
    seu <- seu[keep_gene,]
    print(paste(sample_name,"Filtered gene count:", nrow(seu)))
    return(seu)
  }
}

seurat_preprocess <- function(seu,sct = FALSE,var.to.regress = TRUE,vars =c("percent.mt"),scale.all.gene = FALSE,verbose = TRUE){
  cat("运行前的assay:",DefaultAssay(seu),"\n")
  if (sct == FALSE){
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu,selection.method = "vst", nfeatures = ifelse(2*median(seu$nFeature_RNA) < 5000,2*median(seu$nFeature_RNA),5000))
    if (var.to.regress == TRUE){
      if (scale.all.gene == TRUE){
        seu <- ScaleData(seu,verbose = verbose, vars.to.regress = vars,features =rownames(seu))
      }else{
        seu <- ScaleData(seu,verbose = verbose, vars.to.regress = vars)
      }
    }else{
      if (scale.all.gene == TRUE){
        seu <- ScaleData(seu,features =rownames(seu))
      }else{
        seu <- ScaleData(seu)
      }
    }
    seu <- RunPCA(seu,npcs = 100,verbose = verbose)
    cat("当前的assay:",DefaultAssay(seu),"\n")
    return(seu)
  }else{
    # Normalization
    if(var.to.regress == TRUE){
      seu <- SCTransform(seu, assay = "RNA", variable.features.n = 3000,vars.to.regress = vars) # variable.features.n set to default
    }else{
      seu <- SCTransform(seu, assay = "RNA", variable.features.n = 3000) # variable.features.n set to default
    }
    # Run PCA
    seu <- RunPCA(seu, npcs = 100 ,assay = "SCT", verbose = verbose)
    cat("当前的assay:",DefaultAssay(seu),"\n")
    return(seu)
  }
}

pc_select <- function(seu,reduction = "pca",ndim = 50){
  sample_name <- unique(seu$orig.ident)
  #选择拐点
  a <- ElbowPlot(seu,ndims = ndim,reduction = reduction) + ggtitle(paste0("PC of",sample_name))
  print(a)
  #主成分累积(cumu)贡献>90%
  pct <- seu[[reduction]]@stdev / sum(seu[[reduction]]@stdev) *100
  cat("每个pca主成分的贡献比例如下:","\n",pct,"\n","\n")
  cumu <- cumsum(pct)
  cat("累积主成分的比例如下:","\n",cumu,"\n","\n")
  b <- which(cumu > 90 & pct < 5)[1]
  cat("累积贡献率首次超过90%且单个贡献率小于5%的主成分位置如下:","\n",b,"\n","\n")
  c <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  cat("相邻主成分贡献率差值大于0.1的最后一个pc + 1的位置如下: ","\n",c,"\n","\n")
  pc.select <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1) 
  cat("可能的合适主成分数量为: ",pc.select,"\n")
  return(pc.select)
}

find_neighbors_cluster_umap <- function(seu,reduction = "pca",sct = FALSE,cluster.resns=c(seq(0,1,0.1)),ndim = 30,algorithm = 1,umap.reduction.name ="umap") {
  # original as in using PCA, rather than harmony corrected PCA
  if(sct == FALSE){
    DefaultAssay(seu) <- "RNA"
    seu <- FindNeighbors(seu, 
                         reduction = reduction,dims = 1:ndim) # neighbours are used in FindClusters, but not in RunUmap
    
    seu <- RunUMAP(seu, 
                   reduction = reduction,
                   dims = 1:ndim, 
                   reduction.name = umap.reduction.name)
    
    for (c.res in cluster.resns) {
      seu <- FindClusters(seu, 
                          resolution = c.res,
                          alg=algorithm)
      # seu <- AddMetaData(seu, 
      #                    seu@meta.data$seurat_clusters, 
      #                    col.name = paste0('RNA_snn_', '.', c.res))
    }
    return(seu)
  }else{
    DefaultAssay(seu) <- "SCT"
    seu <- FindNeighbors(seu, 
                         reduction = reduction ,dims = 1:ndim) # neighbours are used in FindClusters, but not in RunUmap
    
    seu <- RunUMAP(seu, 
                   reduction=reduction,
                   dims = 1:ndim, 
                   reduction.name = umap.reduction.name)
    
    for (c.res in cluster.resns) {
      seu <- FindClusters(seu, 
                          resolution = c.res,
                          alg=algorithm)
      
      #seu <- AddMetaData(seu, 
      #                  seu@meta.data$seurat_clusters, 
      #                 col.name = paste0('SCT_snn_', 'res.', c.res))
    }
    return(seu)
  }
}

gsea_df <- read_csv("~/Downloads/lap_vs_dlp_gsea.csv")
gsea_df <- gsea_df[c(1:20,22:27),]
gsea_df <- dplyr::select(gsea_df,ID,enrichmentScore,NES,p.adjust)
gsea_df$reg <- ifelse(gsea_df$NES > 0,"Pos","Neg") 
gsea_df <- gsea_df[order(gsea_df$NES,decreasing = F),]
gsea_df$ID <- factor(gsea_df$ID,levels = gsea_df$ID) #ID因子化排序

p <- ggplot(gsea_df, aes(x = ID, y = NES, fill = reg)) +
  geom_col()  +
  coord_flip() +  # 翻转坐标轴，使条形水平显示
  scale_fill_manual(values = c("blue","red")) +  # 根据正负值填充颜色
  labs(x = NULL, y = "NES", title = "GSEA Enrichment") +
  theme_bw()

#依次从下到上添加标签
up_pathway <- length(which(gsea_df$NES > 0 ))
down_pathway <- length(which(gsea_df$NES < 0 ))
high <- nrow(gsea_df)

p1 <- p + geom_text(data = gsea_df[1:down_pathway,],aes(x = ID,y = 0.1,label = ID),
                    hjust = 0,color = 'black',size = 5) +
  geom_text(data = gsea_df[(down_pathway +1):high,],aes(x = ID,y = -0.1,label = ID),
            hjust = 1,color = 'black',size = 5) +
  scale_x_discrete(labels = NULL) + 
  theme(legend.position = "bottom",
        ##删除网格线与纵坐标轴
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x = element_text(hjust = 0.5,size = 20), 
        axis.ticks.y = element_blank(), ## 删去y轴刻度线
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(vjust = -1, hjust = 0.5, angle = 0, face = "bold", size = 20), 
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1,1,1,1), "cm"),#画布边缘距离上(top)、右(right)、下(bottom)、左(left) 
        plot.title = element_text(hjust = 0.5,size =  22),
        legend.title = element_text(size = 0), 
        legend.text = element_text(size = 22), 
        legend.position = "bottom",
        legend.background = element_rect(fill = 'transparent'))