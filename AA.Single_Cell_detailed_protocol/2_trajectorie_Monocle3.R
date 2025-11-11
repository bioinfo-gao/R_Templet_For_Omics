# https://download.csdn.net/download/weixin_41960569/10864472?spm=1001.2101.3001.6650.6&utm_medium=distribute.pc_relevant.none-task-download-2%7Edefault%7EOPENSEARCH%7ERate-6-10864472-blog-147811138.235%5Ev43%5Epc_blog_bottom_relevance_base6&depth_1-utm_source=distribute.pc_relevant.none-task-download-2%7Edefault%7EOPENSEARCH%7ERate-6-10864472-blog-147811138.235%5Ev43%5Epc_blog_bottom_relevance_base6&utm_relevant_index=12
# https://cole-trapnell-lab.github.io/monocle3/docs/introduction/
# https://cole-trapnell-lab.github.io/monocle3/docs/installation/
# https://zhuanlan.zhihu.com/p/564060683
# https://stuartlab.org/signac/articles/monocle.html

###################  ++++++++++++ =========== ###########################
# 为了避免再次遇到权限警告（针对 C:/Program Files/... 路径），请养成习惯：在安装复杂的、非核心的 $CRAN$ 或 $BiocManager$ 包时，始终使用 lib = "C:/Users/zhen-/Code/ZG_R_Lib" 参数。这能确保 $R$ 只关注您的安全路径，绕过系统库的权限检查。
###################  ++++++++++++ =========== ###########################

# install.packages("devtools") # install.packages("devtools", lib = "C:/Users/zhen-/Code/ZG_R_Lib")
# if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
# 
# devtools::install_github("satijalab/seurat-wrappers")
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'ggrastr'))
# BiocManager::install("Seqinfo", lib = "C:/Users/zhen-/Code/ZG_R_Lib") # BiocManager::install("Seqinfo")
# BiocManager::install("GenomicRanges", lib = "C:/Users/zhen-/Code/ZG_R_Lib", force = TRUE)


# 强制从 GitHub 安装 monocle3，并指定到您的安全路径
# devtools::install_github(
#     "Cole-Trapnell-Lab/monocle3", 
#     lib = "C:/Users/zhen-/Code/ZG_R_Lib"
# )

library("Signac") # library("Signac", lib.loc = "C:/Users/zhen-/Code/ZG_R_Lib")
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)



