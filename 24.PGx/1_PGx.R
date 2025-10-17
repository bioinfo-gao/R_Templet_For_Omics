# 文件：cyp2d6_variant_analysis.R
# 运行环境：RStudio 或 VS Code，Conda r_env（4.3.3）

# 安装包
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install.packages(c("ggplot2", "dplyr"))
setwd("~/Code/R_code/R_Templet_For_Omics/24.PGx")

BiocManager::install("VariantAnnotation")
# 加载库
library(VariantAnnotation)
library(ggplot2)
library(dplyr)

# 步骤 1：模拟 VCF 数据（CYP2D6*4, rs3892097）
vcf_data <- data.frame(
  CHROM = rep("22", 5),
  POS = rep(42128945, 5),  # rs3892097 位置（hg38）
  ID = rep("rs3892097", 5),
  REF = rep("G", 5),
  ALT = c("A", "A", "G", "A", "G"),  # 模拟基因型
  SAMPLE = paste0("Sample", 1:5),
  GT = c("0/1", "1/1", "0/0", "0/1", "0/0"),  # 0=REF, 1=ALT
  DP = c(30, 40, 35, 38, 32),  # 测序深度
  GQ = c(99, 95, 98, 97, 96)   # 基因型质量
)

# 步骤 2：注释 CYP2D6*4 表型
vcf_data <- vcf_data %>%
  mutate(Phenotype = case_when(
    GT == "1/1" ~ "Poor Metabolizer (PM)",  # *4/*4
    GT == "0/1" ~ "Intermediate Metabolizer (IM)",  # *4/NM
    GT == "0/0" ~ "Normal Metabolizer (NM)"  # NM/NM
  ))

# 步骤 3：可视化基因型分布
ggplot(vcf_data, aes(x = SAMPLE, y = Phenotype, fill = Phenotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "CYP2D6*4 表型分布 (rs3892097)",
       x = "样本", y = "表型") +
  scale_fill_manual(values = c("Poor Metabolizer (PM)" = "#EF4444",
                               "Intermediate Metabolizer (IM)" = "#FBBF24",
                               "Normal Metabolizer (NM)" = "#10B981"))

# 保存图表
ggsave("cyp2d6_phenotype_plot.png", width = 8, height = 6)

# 步骤 4：保存 VCF 数据
write.csv(vcf_data, "cyp2d6_vcf.csv", row.names = FALSE)

# 步骤 5：Git 提交
system("git add .")
system("git commit -m 'Add CYP2D6*4 variant analysis'")
system("git push origin main")

# 输出说明
cat("CYP2D6*4 表型分布图保存为 cyp2d6_phenotype_plot.png\n")
cat("VCF 数据保存为 cyp2d6_vcf.csv\n")