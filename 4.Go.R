library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

genes_df <- read.table('genes.txt') 
gene_convert <- bitr(genes_df$V1, fromType = "SYMBOL", 
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org.Hs.eg.db)

markers <- genes_df %>% inner_join(gene_convert,by = c("gene"="SYMBOL"))

all.GO <- enrichGO(gene = deg_up$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)

allgo <- as.data.frame(all.GO)

p <- dotplot(all.GO,label_format = 10000,showCategory = 20)