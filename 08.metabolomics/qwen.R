# 安装必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("xcms", "CAMERA", "pheatmap", "limma", "MSnbase", "ggplot2"))
BiocManager::install("faahKO")




library(MSnbase)
library(xcms)
library(faahKO)
library(MsFeatures)


# 加载R包
library(xcms)
library(CAMERA)
library(pheatmap)
library(limma)
library(ggplot2)

library(MSnbase)
library(xcms)
library(faahKO)
library(MsFeatures)

xmse <- loadXcmsData("xmse")
str(xmse)

# 设置文件路径
file_path <- "path/to/your/raw_data"  # 替换为你的原始数据路径

# 创建xcmsSet对象
raw_data <- xcmsSet(file_path, method = "centWave", ppm = 15, peakwidth = c(5, 20))

# 峰对齐
aligned_data <- group(raw_data, bw = 5, minfrac = 0.5, minsamp = 1)

# 提取峰表
peak_table <- fillPeaks(aligned_data)

# 归一化方法：总峰面积归一化
normalized_data <- sweep(peak_table, 2, colSums(peak_table), "/") * median(colSums(peak_table))

# 或者使用其他归一化方法（如PQN）
source("https://raw.githubusercontent.com/MoseleyBioinformaticsLab/metNormalizer/master/R/metNormalizer.R")
normalized_data <- metNormalizer(peak_table)

# 构建分组信息
group <- factor(c(rep("Control", 5), rep("Treatment", 5)))  # 示例分组

# 使用limma包进行差异分析
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(normalized_data, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取显著差异代谢物
results <- topTable(fit2, adjust = "fdr", number = Inf)
significant_metabolites <- results[results$adj.P.Val < 0.05, ]

# 安装并加载clusterProfiler包
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
library(clusterProfiler)

# 将代谢物ID转换为KEGG ID（假设你有ID映射表）
metabolite_ids <- rownames(significant_metabolites)
kegg_enrichment <- enrichKEGG(gene = metabolite_ids, organism = "hsa", pvalueCutoff = 0.05)

# 可视化通路富集结果
dotplot(kegg_enrichment)

pheatmap(normalized_data, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")

# PCA分析
pca_result <- prcomp(t(normalized_data), scale. = TRUE)
plot(pca_result$x[, 1:2], col = as.numeric(group), pch = 19, xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(group), col = 1:length(levels(group)), pch = 19)

# 保存差异代谢物结果
write.csv(significant_metabolites, "significant_metabolites.csv", row.names = TRUE)

# 保存通路分析结果
write.csv(as.data.frame(kegg_enrichment), "kegg_enrichment_results.csv", row.names = FALSE)