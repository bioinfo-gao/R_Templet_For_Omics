# 安装（首次运行时取消注释）
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GOplot"))

# 加载包
# install.packages("GOplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(dplyr)

#======================= 模拟 300 个差异表达基因（DEGs）
set.seed(123)
deg_data <- data.frame(
  gene = paste0("GENE", 1:300),
  log2FC = c(rnorm(200, mean = 1.5, sd = 0.5), rnorm(100, mean = -1.5, sd = 0.5)),
  pvalue = runif(300, min = 1e-10, max = 0.01)
)


deg_data

# 筛选显著 DEGs（|log2FC| > 1 & p < 0.05）
sig_genes <- deg_data %>% filter(abs(log2FC) > 1 & pvalue < 0.05) %>% pull(gene)

cat("显著差异基因数量:", length(sig_genes), "\n")


# 输出: 显著差异基因数量: 300
# 
# 将基因 Symbol 转换为 Entrez ID（GOplot 要求 Entrez ID）
# 

str(sig_genes)

##  ==================  需要真实数据，以下无法完成
gene_df <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#======================= GO 富集分析（生物过程 BP）
ego <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",          # BP: Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500
)




# 提取 top 10 显著 term
top_terms <- ego@result %>%
  arrange(p.adjust) %>%
  head(10) %>%
  select(ID, Description, GeneRatio, p.adjust)

print(top_terms)

# 准备 GOplot 所需数据格式
# Step 1: 获取每个 term 的基因列表
term2gene <- as.data.frame(ego)
term2gene$geneID <- strsplit(term2gene$geneID, "/")

# Step 2: 展开为长格式
library(tidyr)
term_gene_long <- term2gene %>%
  select(ID, Description, geneID) %>%
  unnest(geneID) %>%
  head(50)  # 限制基因数量避免图形过密

# Step 3: 创建 circ 参数
circ <- circle_dat(
  df = top_terms,           # 富集结果
  gene = term_gene_long     # term-基因对应关系
)

# 绘制圈图
GOChord(
  circ, 
  title = "Top 10 GO Biological Process Enrichment",
  space = 0.001,            # 调整基因间距
  gene.order = "logFC",     # 按 logFC 排序（需提供表达数据）
  gene.size = 3,            # 基因标签大小
  nlfc = 0                  # 不显示 logFC 颜色
)

# 将 log2FC 信息合并到 circ 数据
expr_data <- deg_data %>%
  filter(gene %in% unique(term_gene_long$geneID)) %>%
  select(gene, log2FC) %>%
  rename(geneID = gene)

circ_expr <- merge(circ, expr_data, by = "geneID", all.x = TRUE)

# 绘制带颜色的圈图
GOChord(
  circ_expr,
  title = "GO Enrichment with log2FC Coloring",
  space = 0.001,
  gene.order = "logFC",
  gene.size = 3,
  nlfc = 8,                 # 分 8 段颜色
  lfc.col = c("blue", "white", "red")  # 蓝-白-红渐变
)