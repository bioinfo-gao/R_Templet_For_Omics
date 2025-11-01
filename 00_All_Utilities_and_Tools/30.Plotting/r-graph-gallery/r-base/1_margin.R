# https://r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
# Margins area
# 
#


png("scatter_plot.png", width = 8, height = 6, units = "in", res = 300)

# A demonstration is shown in this dir  
par(oma=c(3,3,3,3))        # for outer margin area. each has 3 lines of space
par(mar=c(5,4,4,2) + 0.1)  # for margin 


#  mai() and omi() set the areas in inches and not in lines.
# 默认值通常为 c(1, 1, 1, 1)。
# 默认值通常为 c(0, 0, 0, 0)。


# Plot
plot(0:10, 0:10, type="n", xlab="X", ylab="Y") # type="n" hides the points


# Place text in the plot and color everything plot-related red
text(5,5, "Plot", col="red", cex=2)
box(col="red")


# Place text in the margins and label the margins, all in forestgreen  
mtext("Margins",             side=3, line=2,          cex=2, col="forestgreen")  
mtext("par(mar=c(b,l,t,r))", side=3, line=1,          cex=1, col="forestgreen")  

mtext("Line 0",              side=3, line=0, adj=1.0, cex=1, col="forestgreen")  
mtext("Line 1",              side=3, line=1, adj=1.0, cex=1, col="forestgreen")  
mtext("Line 2",              side=3, line=2, adj=1.0, cex=1, col="forestgreen")  
mtext("Line 3",              side=3, line=3, adj=1.0, cex=1, col="forestgreen")  

box("figure",   col="forestgreen")  

# Label the outer margin area and color it blue  
# Note the 'outer=TRUE' command moves us from the figure margins to the outer margins.  
# 注意： outer margin 只有保存以后才能看到，在默认图片中不可见
mtext("Outer Margin Area", side=1, line=1, cex=2, col="blue", outer=TRUE)  
mtext("par(oma=c(b,l,t,r))", side=1, line=2, cex=1, col="blue", outer=TRUE)  
mtext("Line 0", side=1, line=0, adj=0.0, cex=1, col="blue", outer=TRUE)  
mtext("Line 1", side=1, line=1, adj=0.0, cex=1, col="blue", outer=TRUE)  
mtext("Line 2", side=1, line=2, adj=0.0, cex=1, col="blue", outer=TRUE)  
box("outer", col="blue")  

# 关闭绘图设备
dev.off()

if(F){ # not needed , demon cm to inch only , see bottom 
  
  # 厘米到英寸转换函数
  cm_to_inch <- function(cm) {
    return(cm / 2.54)
  }
  
  # 设置边距（以厘米为单位）
  margins_cm         <- c(bottom = 2, left = 3, top = 2, right = 1)  # 厘米
  margins_inch       <- cm_to_inch(margins_cm)  # 转换为英寸
  
  outer_margins_cm   <- c(bottom = 1, left = 1, top = 1, right = 1)  # 厘米
  outer_margins_inch <- cm_to_inch(outer_margins_cm)  # 转换为英寸
  
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
}
