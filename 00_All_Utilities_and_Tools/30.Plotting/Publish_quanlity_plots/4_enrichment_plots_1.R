# # 安装必要包（若未安装）
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler",  force = TRUE )
# BiocManager::install("org.Hs.eg.db")
# 
# install.packages(c("ggplot2", "enrichplot"))

# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 步骤 1: 准备示范数据
# 这里使用 clusterProfiler 示例中的基因列表（Entrez ID），模拟一个 GO 富集分析的输入
# 实际课题中替换为您的差异表达基因列表（如来自 DESeq2 或 edgeR 的基因 ID）
gene_list <- c(4312, 8318, 10874, 55143, 79699, 10987, 8407, 6387, 83737, 1234, 5678, 91011)  # 示例 Entrez ID



# 步骤 2: 进行 GO 富集分析
# enrichGO 函数：对基因列表进行 GO 富集（BP: Biological Process 生物过程）
# 参数解释：
# - gene: 输入基因列表 (Entrez ID)
# - OrgDb: 物种数据库 (这里为人 org.Hs.eg.db)
# - keyType: 基因 ID 类型 (Entrez)
# - ont: GO 类型 (BP: 生物过程, CC: 细胞组分, MF: 分子功能)
# - pvalueCutoff: p 值阈值
# - qvalueCutoff: q 值阈值
# - readable: TRUE 将 Entrez ID 转换为基因符号，便于阅读

?enrichGO

data(geneList, package = "DOSE")

de <- names(geneList)[1:100]

go_enrich <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)

# go_enrich <- enrichGO(gene = gene_list,
#                       OrgDb = org.Hs.eg.db,
#                       keyType = "ENTREZID",
#                       ont = "BP",
#                       pvalueCutoff = 0.05,
#                       qvalueCutoff = 0.2,
#                       readable = TRUE)

# 查看富集结果（前 5 行）
go_enrich
head(go_enrich, 5)

# 步骤 3: 绘制 GO 富集圈图 (cnetplot: 概念网络图)
# cnetplot 函数：绘制 GO terms 和基因的网络圈图
# 参数解释：
# - x: enrichGO 结果对象
# - showCategory: 显示前 N 个 GO terms
# - foldChange: 可选，输入基因的 fold change 值（这里模拟）
# - layout: 布局方式 (circle: 圈图)
# - circular: TRUE 启用圈形布局
# - colorEdge: TRUE 根据 p 值着色边缘
foldChange <- rnorm(length(gene_list))  # 模拟 fold change 值

names(foldChange) <- go_enrich$geneID[1:length(gene_list)]  # 匹配基因 ID


cnetplot(go_enrich, showCategory = 10, foldChange = foldChange, layout = "circle", circular = TRUE, colorEdge = TRUE)

# 保存圈图到文件
ggsave("go_enrich_circle_plot.png", width = 10, height = 10)

# 步骤 4: Git 提交（需在终端运行或通过 system 命令）
# 先确保稀疏检出包含路径
# system("git sparse-checkout add 'Finished Versions/Ch01/01_03'")
# system("git add go_enrich_circle_plot.png")
# system("git commit -m 'Add GO enrichment circle plot'")
# system("git push origin main")

# 输出说明
cat("GO 富集圈图已绘制并保存为 go_enrich_circle_plot.png\n")
cat("示例基因列表: ", paste(gene_list, collapse = ", "), "\n")
cat("GO 富集结果: ", nrow(go_enrich), " 个 terms\n")
