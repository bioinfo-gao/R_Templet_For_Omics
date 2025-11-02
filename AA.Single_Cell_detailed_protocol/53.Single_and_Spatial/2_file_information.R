# 文件：scRNA_seq_reads.R
# 运行环境：RStudio 或 VS Code，Conda r_env（4.3.3）
setwd("~/Code/R_code/R_Templet_For_Omics/23.Single_and_Spatial")
# 加载库
library(Seurat)
library(ggplot2)

# 模拟读取I_1和R_2文件（实际需替换为真实FASTQ路径）
# 假设已运行Cell Ranger，生成表达矩阵
# PBMC3k 数据路径：http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")

# 创建Seurat对象
pbmc_seurat <- CreateSeuratObject(counts = pbmc, project = "PBMC3K", min.cells = 3, min.features = 200)

# 模拟I_1样本索引（示例数据）
sample_index <- data.frame(
  Sample = c("Sample1", "Sample2"),
  Index_I1 = c("ATCGTAGC", "GCTACCTA")
)

# 模拟R_2基因表达统计
gene_counts <- Matrix::rowSums(pbmc_seurat@assays$RNA@counts)
top_genes <- names(sort(gene_counts, decreasing = TRUE))[1:5]

# 可视化R_2基因表达
pdf("r2_gene_counts.pdf")
barplot(gene_counts[top_genes], names.arg = top_genes, col = "blue",
        main = "Top 5 Genes from R_2 (PBMC3K)", xlab = "Genes", ylab = "Counts")
dev.off()

# 保存I_1索引数据
write.csv(sample_index, "i1_index.csv", row.names = FALSE)

# Git 提交（需在终端运行）
# git sparse-checkout add "Finished Versions/Ch01/01_03"
# git add "Finished Versions/Ch01/01_03/i1_index.csv"
# git add "Finished Versions/Ch01/01_03/r2_gene_counts.pdf"
# git commit -m "Add I_1 and R_2 analysis for scRNA-seq"
# git push origin main

# 输出说明
cat("I_1索引保存为 i1_index.csv\n")
cat("R_2基因表达图保存为 r2_gene_counts.pdf\n")