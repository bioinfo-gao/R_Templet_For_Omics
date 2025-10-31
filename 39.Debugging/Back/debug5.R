# 文件：debug_cyp2d6_phenotype.R
# 运行环境：RStudio，Conda R43 (R 4.3.3)

# 步骤 1：定义函数（模拟 CYP2D6*4 表型计算）
calculate_phenotype <- function(gt) {
    if (gt == "1/1") {
        phenotype <- "Poor Metabolizer (PM)"  # 调试断点：检查 gt 值
    } else if (gt == "0/1") {
        phenotype <- "Intermediate Metabolizer (IM)"  # 调试断点：逐步执行
    } else {
        phenotype <- "Normal Metabolizer (NM)"
    }
    return(phenotype)
}

# 步骤 2：测试数据
gt_list <- c("0/1", "1/1", "0/0", "0/1", "0/0")

# 步骤 3：调用函数
phenotypes <- sapply(gt_list, calculate_phenotype)

# 步骤 4：可视化
library(ggplot2)
df <- data.frame(GT = gt_list, Phenotype = phenotypes)
ggplot(df, aes(x = GT, y = Phenotype, fill = Phenotype)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "CYP2D6*4 表型分布 (rs3892097)")

# 保存图表
ggsave("cyp2d6_debug_plot.png", width = 8, height = 6)


# 输出
cat("图表保存为 cyp2d6_debug_plot.png\n")