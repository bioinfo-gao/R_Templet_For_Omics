# 创建示例字符串向量（模拟 CAH 基因数据）
genes <- c("CYP2C19*1", "CYP2C19*2", "CYP2C19*1", "CYP2C19*3", "CYP2C19*2")

# 统计频数
freq_table <- table(genes)
print(freq_table)

library(tidyverse)

# 创建数据框
data <- tibble(genes = c("CYP2C19*1", "CYP2C19*2", "CYP2C19*1", "CYP2C19*3", "CYP2C19*2"))

# 统计频数
freq_count <- data %>% count(genes)
print(freq_count)