# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# 确保 BiocManager 已安装
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 强制重新安装所有核心 Bioconductor 依赖
BiocManager::install(c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges"), 
                     ask = FALSE, 
                     force = TRUE)

# 1. 定义一个新的、干净的临时目录路径（例如，直接在 C 盘根目录下）
# 确保这个目录事先不存在
new_temp_dir <- "C:/R_Temp_Install" 

# 2. 创建该目录
if (!dir.exists(new_temp_dir)) {
    dir.create(new_temp_dir, recursive = TRUE)
}

# 3. 设置 R 的环境变量，将所有临时文件和下载路径指向这个新目录
Sys.setenv(TMPDIR = new_temp_dir)
Sys.setenv(TEMP = new_temp_dir)
Sys.setenv(R_INSTALL_STAGED = FALSE) # 确保直接写入，避免复杂路径

message(paste("R 的临时安装路径已设置为:", new_temp_dir))

BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86"), 
                     ask = FALSE, 
                     update = FALSE)



# 重新安装 DESeq2，使用 type="binary" 避免源代码编译错误
BiocManager::install("DESeq2", ask = FALSE, type = "binary")


install.packages("devtools")
devtools::install_github("satijalab/azimuth")

# 验证核心依赖是否修复
library(GenomeInfoDb)
library(GenomicRanges)
library("Azimuth")

# Install the PBMC systematic comparative analyis (pmbcsca) dataset
InstallData("pbmcsca")

# returns a Seurat object named pbmcsca
pbmcsca <- LoadData("pbmcsca")

# The RunAzimuth function can take a Seurat object as input
pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")


# We can now visualize the outputs of Azimuth. Note that all cells are labeled with high-resolution annotations, and are projected into a harmonized space despite being analyzed from a wide variety of technologies.

p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmcsca, group.by = "Method")
p1 + p2