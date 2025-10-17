# 文件：plot_margin_cm.R
# 运行环境：RStudio 或 VS Code，Conda r_env（4.3.3）

# 加载库
library(Seurat)
library(ggplot2)

# 厘米到英寸转换函数
cm_to_inch <- function(cm) {
  return(cm / 2.54)
}

# 设置边距（以厘米为单位）
margins_cm         <- c(bottom = 2, left = 3, top = 2, right = 1)  # 厘米
margins_inch       <- cm_to_inch(margins_cm)  # 转换为英寸
outer_margins_cm   <- c(bottom = 1, left = 1, top = 1, right = 1)  # 厘米
outer_margins_inch <- cm_to_inch(outer_margins_cm)  # 转换为英寸



#install.packages('devtools')
devtools::install_github('satijalab/seurat-data')
library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)

InstallData("pbmc3k")  
data("pbmc3k")  

pbmc <-  pbmc3k.final # sce

# 加载 PBMC3k 数据（Seurat 示例数据）
#pbmc <- readRDS(url("http://satijalab.org/seurat/pbmc3k_final.rds"))  # 需联网

# 提取基因表达数据（示例：MS4A1 和 CD3D）
data <- as.data.frame(t(as.matrix(pbmc@assays$RNA@data[c("MS4A1", "CD3D"), ])))
colnames(data) <- c("MS4A1", "CD3D")

# 设置绘图设备
png("scatter_plot_cm_margins.png", width = 8, height = 6, units = "in", res = 300)

# 设置边距
par(mai = margins_inch, omi = outer_margins_inch)

# 绘制散点图
plot(data$MS4A1, data$CD3D,
     xlab = "MS4A1 Expression", ylab = "CD3D Expression",
     main = "Gene Expression Scatter Plot with Custom Margins (cm)",
     pch = 16, col = "blue", cex = 0.5)

# 关闭绘图设备
dev.off()

getwd() 
setwd("C:/Users/zhen-/Documents/Code/Rcode/R_Templet_For_Omics/plotting/r-graph-gallery/r-base/") 

# Git 提交（需在终端运行）
# git sparse-checkout add "Finished Versions/Ch01/01_03"
# git add "Finished Versions/Ch01/01_03/scatter_plot_cm_margins.png"
# git commit -m "Add scatter plot with cm margins"
# git push origin main

# 输出说明
cat("散点图已保存为 scatter_plot_cm_margins.png\n")
cat("内边距 (cm): ", paste(margins_cm, collapse = ", "), "\n")
cat("外边距 (cm): ", paste(outer_margins_cm, collapse = ", "), "\n")
