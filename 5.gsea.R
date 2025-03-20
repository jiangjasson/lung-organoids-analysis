library(clusterProfiler)
library(msigdf)
library(enrichplot)

Idents(sce_filter) <- sce_filter$celltype
sce_filter <- JoinLayers(sce_filter)
sce_filter_all_marker <- FindMarkers(sce_filter,only.pos = FALSE,
                                     ident.1 = "LAP",ident.2 = "d-LPC",
                                     assay = 'RNA',slot = 'data',
                                     min.pct = 0,
                                     logfc.threshold = 0) %>% mutate(difference = .$pct.1 - .$pct.2, gene = rownames(.)) 
#更具fc排序
genelist <- sce_filter_all_marker$avg_log2FC
names(genelist) <- rownames(sce_filter_all_marker) 
genelist <- sort(genelist,decreasing = T)

a <- msigdbr(species = "Homo sapiens",category = "C5")
a <- a %>% filter(gs_subcat != "HPO") %>% dplyr::select(gs_name,gene_symbol)

msig_gsea <- clusterProfiler::GSEA(genelist,  #msiDB数据库输入的是symbol
                                   TERM2GENE = a,
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   pvalueCutoff = 0.05,
                                   nPermSimple = 100000,
                                   eps = 0
)

msig_gsea_result <- msig_gsea@result
write.csv(msig_gsea_result,file = "OUT/gsea_go.csv")

gsea_df <- read_csv("gsea.csv")
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