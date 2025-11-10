# 安装 Seurat (如果还没有)
# install.packages("Seurat")

# 安装其他依赖包
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'GO.db', 'limma'))

# 加载 Seurat 和其他绘图包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/BB_Spatial_Transcriptomics")
# 下载示例数据 (如果需要，这是一个大约 50MB 的文件)
# download.file("https://cf.10xgenomics.com/samples/spatial-gene-expression/1.0.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz", destfile = "V1_Adult_Mouse_Brain.tar.gz")
# untar("V1_Adult_Mouse_Brain.tar.gz")

# 使用 Load10X_Spatial 函数加载数据
# "filtered_feature_bc_matrix" 是 Space Ranger 的标准输出目录名
slice <- Load10X_Spatial(data.dir = "C:/Users/zhen-/Code/R_code/R_Omics_DS/BB_Spatial_Transcriptomics/filtered_feature_bc_matrix/")

# 查看 Seurat 对象
print(slice)
# An object of class Seurat
# ... features, ... cells
# 1 spatial image assay
