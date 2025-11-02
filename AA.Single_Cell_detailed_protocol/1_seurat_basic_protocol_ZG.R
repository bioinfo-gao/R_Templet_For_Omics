#  https://rpubs.com/shengwei_lin/1082110
# https://blog.csdn.net/weixin_43843918/article/details/141168265  Best 
# https://zhuanlan.zhihu.com/p/14221122458
# https://cloud.tencent.com/developer/article/2224200
# https://blog.csdn.net/2301_78336013/article/details/147811138 
rm(list = ls()) 
setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/AA.Single_Cell_detailed_protocol")
library(dplyr)
library(Seurat)
library(patchwork)

################## ============  01 Read Data  ============ ################## 
# pbmc : Peripheral blood mononuclear cells
pbmc.data <- Read10X(data.dir = "./53.Single_and_Spatial/filtered_gene_bc_matrices/hg19/") # pbmc

# 创建一个Seurat对象:
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells =3 , min.features = 3)

#查看对象信息
pbmc # 32738 genes ==> 13714, cells 2700 ==> 2700 
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


######### ######### ============  02 标准预处理工作流程  ====== ====== ##################
# Seurat中一些常用的质控(QC)指标:
# 
# nFeature_RNA:每个细胞中检测到的基因的数量。 # 死细胞，细胞碎片或空液滴通常很少基因; 双/多胞体可表现出异常高的基因数。
# nCount_RNA:细胞内检测到的分子总数。         # UMI的总数量(一个基因可能会因为有多个转录本而被检测到多次),
# percent.mt:对应线粒体基因组的序列的百分比。 # 低质量/濒死的细胞通常表现出广泛的线粒体污染。

head(pbmc@meta.data)
row.names(pbmc)[grep("^MT-", rownames(pbmc)) ] #rownames = row.names   Total : 13  gene 

# "[[]]"运算符可以向对象pbmc的元数据meta.data中添加列。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")  
# 对于未整理好的线粒体基因，直接安装线粒体基因的通用名，总共只有13个基因
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "ATP6|ATP8|COX1|COX2|COX3|CYTB|ND1|ND2|ND3|ND4|ND4L|ND5|ND6")

VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA",  "percent.mt"),  ncol = 3)

# featureScatter是一种典型的可视化特征与特征之间关系的方法
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

# 根据小提琴图和散点图的情况，个人认为可以选择如下的质控过滤指标，
# 即： nFeature_RNA > 200（最开始的min.features = 200）& nFeature_RNA < 2000 nCount_RNA < 7000 percent.mt < 7 
# 但是呢，为了和官方教程相统一，我们还是采用官方的过滤阈值，实际上相差并不多! 官方 {200 2500} 5
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200  & 
                   nFeature_RNA < 2500 & 
                   nCount_RNA   < 5000 &
                   percent.mt   < 5 )
pbmc 

# check data after QC
VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA",  "percent.mt"),  ncol = 3)
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt"  )
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3+plot4


################## ====== ======  03 数据标准化  ====== ====== ##################
## 数据标准化 通过质控过滤完不需要的细胞之后呢，就要对数据进行标准化处理了！
#  default： LogNormalize：全局缩放的标准化方法（其实还有另一种标准化方法SCTransform） 
#  数据标准化的目的便是排除技术误差，让测序深度和文库大小不同的细胞的基因表达量具有可比性，
# log1p指的是log(x + 1)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
head(GetAssayData(pbmc, assay = "RNA", slot = "data")[1:6, 1:25])
# 在 Seurat v5 (Assay5 对象) 中，内部数据的存储方式发生了变化，
# 如果您想访问原始计数 (Raw Counts)
head(GetAssayData(pbmc, assay = "RNA", slot = "counts")[1:6, 1:25])


################## ====== ======  04 鉴定高可变基因 ================ ################## 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) 
# 鉴定10个最高变的基因
top10 <- head(VariableFeatures(pbmc), 10)
# 绘制高可变基因的火山图
plot5 <- VariableFeaturePlot(pbmc)
plot6 <- LabelPoints(plot= plot5, points= top10, repel = TRUE) # ggrepel too label gene
plot5 + plot6 # Visible Only when you make the plot area very huge

combined_plot <- plot5 + plot6
# pdf(file=paste0(gene, ".PFS.pdf"), width=5.5, height=5, onefile=FALSE)
pdf("expression_and_variance.pdf", width = 12, height = 6)
combined_plot# << ========== plot is huge NOT visible in regular window
dev.off()


################## ============  05 数据缩放（数据中心化）   ============ ##################  
# 调整每个基因的表达量，使得细胞间的平均表达量为0 缩放每个基因的表达量，使得细胞间的方差为1 目的：这一步骤赋予每个基因在下游分析中同样的权重，
# 不至于使高表达基因占据主导地位 结果存储在pbmc[[“RNA”]]$scale.data中

all.genes <- rownames(pbmc)
# 以基因为单位对矩阵数据进行缩放
pbmc <- ScaleData(pbmc, features = all.genes) # features= scale.genes  will be faster  



################## ============  06 PCA   ============ ##################  
pbmc <- RunPCA(pbmc, features = VariableFeatures (pbmc))

# PCA结果的可视化（其实不太重要）。 
# Seurat提供了几种方法来对PCA结果的细胞和基因进行可视化，包括VizDimReduce()、 DimPlot()和DimHeatmap()。 

# 第1种：VizDimLoadings()：它显示了每个变量在每个主成分上的贡献程度，
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #没啥用其实

# 第2种：DimPlot()：（PCA, t-SNE, UMAP） ==>> 第一个线性降维，第二第三非线性降维
DimPlot(pbmc, reduction = "pca") # #也没啥用,

# 第3种：DimHeatmap()
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # 也没啥用
DimHeatmap(pbmc, dims = 1,    cells = 500, balanced = TRUE) # 也没啥用


################## ============  07 确定后续分析需要的“维度”）   ============ ##################  

ElbowPlot(pbmc) # 一般选择肘点附近（即本图中PC=10左右），作为数据集的PC值

# ElbowPlot from JackStraw  =============》》 不是必需，想看看一看
# JackStraw 分析是一种用于评估主成分分析中每个主成分显著性的方法，通过随机重采样的方式来模拟随机背景分布，
# 然后将真实的主成分与随机背景进行比较，从而判断每个主成分是否包含真实的生物学信号。 

# num.replicate = 100：指定随机重采样的次数，
# pbmc=JackStraw(pbmc, num.replicate = 100) 
# pbmc=ScoreJackStraw(pbmc, dims = 1:15)
# JackStrawPlot(pbmc,dim = 1:5)

# 横坐标的PC是你接下来要选择分析的维度数，即多少组基因表达pattern。当选择一组的时候，即所有基因都在一个表达pattern里面（我们从生物学角度看也肯定是不可能的），我们很容易想到，
# 所有细胞都在一个组内，肯定是存在很大差异的，所以方差值（方差=标准差（Standard Deviation）的平方）就会很大。当PC增多，即分组增多后，
# 每个组内的差异也就慢慢变小了，方差值自然也就变小了。再往后，随着选的PC的增多，基本上区分不出来组内的差异了，这个时候方差基本就稳定了下来，
# 此时在选取更大的PC值也就没有更大的意义了，并且还会增加计算压力。
# 因此，个人认为关于PC选取的最合适原则是在拐点右侧3-5个维度即可（本数据拐点在7-8左右，我们选10刚刚好），宽松一点的话，可以设置在拐点右侧10以内。
#不过最重要的一点还是尽可能地多尝试，选取最适合自己的PC值，这一点也是官网教程所建议的。



################## ============  08 Cell Clustering   ============ ##################   08 细胞聚类
# Seurat 应用了一种基于图表的聚类方法。重要的是,驱动聚类分析的距离度量(基于先前识别的PC) 保持不变。
# 将细胞产制边界,然后将该生均中如K最近邻(KNN)图,在具有相似特征表达模式的细胞之间绘美的不同簇。
# 首先基于PCA空间中的欧几里得度量构建一个KNN图,然后基于局部邻域中的共享重叠(Jaccard 相似性)来精确任意两个细胞之间的边缘权重。
# 这个步骤使用 FindNeighbors()函数执行,并将先前定义的前10个PC作为输入。

pbmc <- FindNeighbors(pbmc, dims = 1:10) 

#为了将细胞聚类，我们接下来应用模块化优化技术，如Louvain算法(默认值)或SLM，迭代地将细胞聚类在一起。
pbmc <- FindClusters(pbmc, resolution=0.5) # 解析参数resolution，可以理解为就是分辨率 value setting: 0.4-1.2



################## ============  09 非线性降维（e.g. UMAP/tSNE）   ============ ##################   08 细胞聚类

# 关于UMAP和TSNE图的展示简单做一个描述：可以说UMAP展现的距离更真实，TSNE图展示的结果似乎更美观
pbmc<-RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap") # 0-7 , total 8 cluster

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")


################## ============  10 寻找表达差异基因   ============ ##################   08 细胞聚类
# 计算所有簇的markers Seurat教程仅设置了一个参数，就是only.pos，该参数的意思是只保留正标记（positive markers），也就是只保留相对于其他簇而言，目标簇显著高表达的基因标记。 
# FindAllMarkers() 函数：这是 Seurat 包中的一个函数，用于对 Seurat 对象 pbmc 中的所有细胞聚类进行差异表达分析，找出每个聚类的标记基因。
# only.pos = TRUE：该参数指定只返回在特定聚类中高表达的基因，即正差异表达的基因，忽略那些在其他聚类中高表达的基因。 
# -min.pct = 0.25：是一个过滤条件，要求目标基因至少在 25% 的当前聚类细胞或至少 25% 的其他聚类细胞中表达，
# 才会被纳入差异表达分析，这样可以减少在极少数细胞中表达的基因对分析结果的干扰。 logfc.threshold = 0.25：
# 表示只考虑平均对数 2 倍变化值大于 0.25 的基因，即该基因在当前聚类中的表达量至少是其他聚类的 2^0.25 = 1.19 倍，这有助于筛选出表达差异更显著的基因。


# 找出每个细胞簇的标记物,与所有剩余的细胞进行比较,只报告阳性细胞
#pbmc.markers<-FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.25)
pbmc.markers<-FindAllMarkers(pbmc, only.pos = TRUE) # 10 min

# just group 
pbmc.markers %>%
    group_by(cluster)

head(pbmc.markers %>% group_by(cluster))
tail(pbmc.markers %>% group_by(cluster))

# each group show the 
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

#slice_max(n = 2, order_by = avg_log2FC)：对于每个分组（即每个聚类），按avg_log2FC 列的值进行降序排序，并选取前 2 行，也就是每个聚类中 avg_log2FC 最高的前 2 个基因。


# 寻找cluster2 的相对于所有其他簇的 markers
cluster2.markers <- FindMarkers(pbmc, ident.1 =2, min.pct = 0.25)
#寻找cluster5 与cluster0和cluster3 不同的所有markers
cluster5.markers <- FindMarkers(pbmc, ident.1= 5, ident.2 = c(0, 3) , min.pct = 0.25)


# 实际上，到了这里单细胞Seurat流程的分析基本上走完了，在这里保存RDS是比较合适的，
# 因为里面也存好了计算的差异基因的情况，到时候重新加载RDS文件时，不需要运算即可获取。
#同时与注释获得完整的细胞类型还有一定的时间要求。所以，我们选择在这里保存第一版分析数据的RDS。
#保存结果为rds文件


output_dir <- "./output"      # dirname(output_file) output_file <- "./output/pbmc_tutorial.rds"

# 检查文件夹是否存在，如果不存在则创建
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("已创建文件夹:", output_dir, "\n")
} 

#  保存 Seurat 对象
saveRDS(pbmc, file = file.path( output_dir, "pbmc_tutorial.rds" ) ) # another method paste0()
cat("Seurat 对象已成功保存到:", output_file, "\n")


################## ============  11  基因表达可视化  ============ ################## 

# 第1种：VizDimLoadings()
# 第2种：DimPlot()：（PCA, t-SNE, UMAP）
# 第3种：DimHeatmap() 
# 以上三种已经在数据缩放 == >> PCA 中使用过了

# 第4种：VlnPlot()   这里继续深入
# 第5种：山脊图      # 第5种实际上是第4种旋转90度的图像
# 第6种：featurePlot
# 第7种：dotplot

# 按簇展示特定基因表达水平
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) 
# 行数据 ==> slot == counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# Above ==>>　X 取离散整数值时, log(X+1) 也会取固定的、离散的值：
# UMI 原始计数 XLog 转换值 log(X+1) (以 e 为底) 
# log(1)  大量的细胞聚集在 Y=0 的线上 
# log(2) 所有 UMI=1 的细胞聚集在 Y≈0.693 的线上
# log(3)1.099$所有 UMI=2 的细胞聚集在 Y≈1.099 的线上
# log(4)1.386$所有 UMI=3 的细胞聚集在 Y≈1.386 的线上
# log(11)2.398$所有 UMI=10 的细胞聚集在 Y≈2.398 的线上
# log(101)4.615$所有 UMI=100 的细胞聚集在 Y≈4.615 的线上


features <- c("LYZ", "CCL5", "IL32", "PTPRCAP")
RidgePlot(pbmc, features = features, ncol = 2)
#Dicking joint handwidth of A 319

# 我命名其为表达特征定点图
FeaturePlot(pbmc, features = features) 


# 点阵图 == 简单
DotPlot(pbmc, features = features) + RotatedAxis()

# 点阵图 == 高级
DotPlot(pbmc, features = list(
    "T-like NKT" =
        c("CD3E","TRAC", "CD8A", "IL7R","KLF2"),
    "NK-like NKT" =
        c("NKG7", "GNLY", "PRF1","KLRD1","CCL5"),
    "Activated" = c("FOS", "JUN", "CD69"),
    "Exhausted" = c("PDCD1","LAG3","TIGIT")
)) + RotatedAxis()

DoHeatmap(subset(pbmc, downsample =100), features = features,  size = 3)

#DoHeatmap()为给定的细胞和特征生成一个表达热图。在这种情况下,我们为每个集群绘制前20个标记(或所有标记,如果小于20)。
pbmc.markers%>%
    group_by(cluster)%>%
    top_n( n = 5, wt=avg_log2FC) -> top5

top5
DoHeatmap(pbmc, features = top5$gene) + NoLegend()

#这里记住一点，基因表达量的图是用来显示基因在各簇中的表达情况的，这个为该基因是否可以作为本簇的marker提供了直接的证据来源，
# 如和已有报道相互印证（如XXX细胞特异性表达该基因），则可通过该基因对某一簇细胞进行细胞类型注释。细胞类型注释是一个人工主观的过程（自动注释的应用范围很少，只聚焦在极少数物种中，
# 大多数的细胞类型注释还是依赖人工进行的），这一过程已游离于上述前10步分析之外，
# 所以这就是为什么我在第10步结束的时候保存.rds的原因，不过这里也确实没有细讲如何进行细胞类型注释，咱们后续专门出一期和大家唠唠细胞注释的内容！ 

################## ============  12 细胞簇命  ============ ################## 
# 假定我们走完上述11步，已经完成了细胞类型注释工作，即知道了每个簇对应的细胞类型。那么现在我们要如何将这些簇更改成对应的细胞名称呢？看我手法：
# 如果您的聚类 ID 数量很多，建议使用 data.frame 或 Excel 表格来整理映射关系，然后导入 R 中转换为命名向量。
#把每个统对应的新名字写成一个能合
#这里注意,如果有两个候注释到的名字一样,
#要在对应数字位置都写一遍<<==即一定保证这个集合的元素数和保的个数是一致的8
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Мемогу CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC") #, "platelet")

new.cluster.ids

#这一步是把pbmc的levels传输给new.cluster.ids作为他的name属性
# 1. 创建命名向量
# 格式为：New Name = Old ID 
# 注意：在 R 中创建命名向量时，通常是将 '值' 赋给 '名称'

names(new.cluster.ids) <- c("0", "1", "2", "3", "4", "5", "6", "7") 
#names(new.cluster.ids) <- c("1", "2", "3", "4", "5", "6", "7", "8","9") 

# 2. 执行重命名操作
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 3. 验证结果
head(Idents(pbmc)) # 查看前几个细胞的新 ID


3 #最后,再展示一下注释后的细胞图谱
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size= 0.5) + NoLegend()


# 到这里，Seurat的基础分析算全部完成了，这时候我们可以保存一下最完整的信息。
saveRDS(pbmc, file="../output/pbmc3k_final.rds")
