library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(clustree)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggvenn)
set.seed(123)

#1.Cell type proportion
metadata <- sce_filter@meta.data
##调整细胞群的位置
metadata$celltype <- factor(metadata$celltype,levels = rev(c("D-LPC","p-LPC","Pre p-LPC","Hindgut","Neural-like","Neuroendocrine","Cycling","Unknown","AFE")))
color_proportion <- c(
  "#33a02c", # 绿色
  "#cab2d6", # 淡紫色
  "#fdbf6f", # 淡橙色
  "#a6cee3", # 淡蓝色
  "#fb9a99", # 淡红色
  "#b2df8a", # 淡绿色
  "#1f78b4", # 蓝色
  "#ff7f00", # 橙色
  "#e31a1c"# 红色
)
# 计算每个样本中各个群体的细胞数目和比例
sample_col <- sym("orig.ident")
group_col <- sym("celltype")
cell_counts <- metadata %>%
  group_by(!!sample_col, !!group_col) %>%
  summarise(count = n()) %>%
  group_by(!!sample_col) %>%
  mutate(proportion = count / sum(count))

p <- ggplot(data = cell_counts, aes(x = !!sample_col, y = proportion, fill = as.factor(!!group_col))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 8),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(size = 16))+
  scale_fill_manual(values = color_proportion) +
  labs(x = NULL, y = "Cell type proportion", fill = "Celltype") 

#2.volcano_plot
df <- read.csv("deg.csv")
volcano_plot_enhanced(df,logFC = "avg_log2FC",FDR = "p_val_adj",Symbol = "gene",gene_colum_name = "gene")

#3.dimplot
main_cell_type_color <- c('#AC142E','#DE901C','#0D5EA4') #,"#64cccf")
p1 <- dimplot_cellnum(sce_filter,group.by = "main_celltype",reduction = "harmony.umap",cols = main_cell_type_color,pt.size = 0.5 )+
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 20),        # 调整图例文本大小
    legend.title = element_text(size = 20),       # 调整图例标题大小
    legend.key.size = unit(1, "cm"),              # 调整图例键的大小
    legend.position = "right"                     # 调整图例的位置
  )
ggsave(filename = "main_celltype.pdf",p1,width = 14,height = 12)

main_cell_type_gene <- c("EPCAM","CDH1","KRT18","COL6A2","COL3A1","PRRX1","DCX","NES","MAP2")
p2 <-  DotPlot(sce_filter,group.by = "main_celltype",features = main_cell_type_gene,cols = 'Spectral' )+ 
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        axis.text.y = element_text(margin = margin(r = 4))) +
  labs(x = NULL, y = NULL)
ggsave(filename = "main_celltype_markergene.pdf",p2,width = 7,height = 4)

#4.vlinplot
sce_temp$orig.ident <- factor(sce_temp$orig.ident,levels = c("HAWO_3D_D35","HAWO_SC_D35","HAWO_SC_D90","HALO_SC_D35","HALO_SC_D90"))
col <- c("#94c58f","#FF9D9A", "#ff0000","#9cd2ed","#0D63A5")
p <- VlnPlot(subset(sce_temp,main_celltype == "LAP"),features = c("NKX2-1","SCGB3A2","SFTPB","CFTR","SLC4A4","MUC1","STEAP4","CEACAM6"),
             group.by = "orig.ident",pt.size = 0,,cols = col,ncol = 2) * 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 8)) * 
  theme_test() * labs(x = NULL) 

#5.heatmap
expression_data <- GetAssayData(sce_filter, layer = "data")[deg_top,]
metadata <- sce_filter@meta.data

metadata <- metadata[order(metadata$pseudotime), ]
# 根据排序后的元数据对表达矩阵重新排序
ordered_expression_data <- expression_data[, rownames(metadata)]

# 定义注释
annotationcol <- data.frame(CellType = metadata$custom_celltype,Pseudotime = metadata$pseudotime)
rownames(annotationcol) <- rownames(metadata)

#为 Pseudotime 创建颜色向量
# 创建颜色渐变函数
pseudotime_colors <- colorRamp2(
  breaks = c(min(annotationcol$Pseudotime), median(annotationcol$Pseudotime), max(annotationcol$Pseudotime)),
  colors = c('#0B1BE9', '#F19F04', 'yellow')
)

# 定义颜色
annotation_colors <- list(
  CellType = c("d-LPC" = '#007ABA', "LAP" = '#E91E25', "AT2" = '#DF75AE'),
  Pseudotime = pseudotime_colors
  #Pseudotime = colorRampPalette(c('#0B1BE9','#F19F04', "yellow"))(100)  # 假设 Pseudotime 是一个连续变量
)

# 定义列注释
ha_col <- HeatmapAnnotation(
  CellType = annotationcol$CellType,
  Pseudotime = annotationcol$Pseudotime,
  col = annotation_colors,
  annotation_legend_param = list(
    Pseudotime = list(at = c(min(annotationcol$Pseudotime), ceiling(median(annotationcol$Pseudotime)),ceiling(max(annotationcol$Pseudotime))))
  )
)

p <- Heatmap(t(scale(t(ordered_expression_data))),
              name = "expression",
              col = colorRamp2(c(-2, 0, 2), c('skyblue1', "grey10", "yellow")),
              cluster_rows = TRUE,
              show_row_dend = FALSE,
              row_dend_reorder = FALSE,
              row_names_side = "left",
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 6),
              top_annotation = ha_col,
              use_raster = FALSE
)


#6.venn plot
venn_list <- list(
  "rsv_Low_LAP_vs_con_Low_LAP" = rsv_low_LAP_vs_con_Low_LAP$gene,
  "High_rsv_LAP_vs_Low_rsv_LAP" = High_rsv_LAP_vs_Low_rsv_LAP$gene,
  "High_rsv_LAP_vs_rsv_non_epi" = High_rsv_LAP_vs_rsv_non_epi$gene,
  "High_rsv_LAP_vs_rsv_PNEC" = High_rsv_LAP_vs_rsv_PNEC$gene
)

p <- ggvenn(venn_list,show_elements = F,
            fill_color = c('#F37121','#0D63A5',"#b2db87","#7ee7bb"),
            fill_alpha = 0.8,
            stroke_color = "#eeeeee",stroke_size = 0,
            set_name_color = "#161515",set_name_size = 4,
            text_color = "#3d293f",text_size = 4)

#7.rsv featureplot
my_colors <- colorRampPalette(c("#00008B","cyan",  "lightgreen","yellow", "orange","red"))(100)

sce_filter$percent.rsv <- log10(sce_filter$percent.rsv+1)

p <- FeaturePlot(object = sce_filter, features = "percent.rsv", pt.size = 1, order = T) + 
  scale_colour_gradientn(
    colors = my_colors
  )

df_rsv <- subset(sce_filter,orig.ident == "d82-rsv")
df_con <- subset(sce_filter,orig.ident == "d82-con")
p1 <- FeaturePlot(object = df_rsv, features = "percent.rsv", pt.size = 1, order = T) + 
  scale_colour_gradientn(
    colors = my_colors
  )
p2 <- FeaturePlot(object = df_con, features = "percent.rsv", pt.size = 1, order = T) + 
  scale_colour_gradientn(
    colors = my_colors
  )