#  https://blog.csdn.net/weixin_43843918/article/details/141168265  Best 
# https://blog.csdn.net/2301_78336013/article/details/147811138 
setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/AA.Single_Cell_detailed_protocol")
library(dplyr)
library(Seurat)
library(patchwork)

rm(list = ls()) 
# pbmc : Peripheral blood mononuclear cells
pbmc.data <- Read10X(data.dir = "./53.Single_and_Spatial/filtered_gene_bc_matrices/hg19/") # pbmc

# 接下来,我们就来创建一个Seurat对象:

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells =3 , min.features = 3)

#查看对象信息
#
pbmc
# 32738 genes ==> 13714, cells 2700 ==> 2700 
# 
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 密集矩阵
dense.size <- object.size(as.matrix(pbmc.data))
dense.size #709591472 bytes
#稀疏矩阵
sparse.size <- object.size(pbmc.data) 
sparse.size #29905192 bytes
# 两者相差倍数
dense.size/sparse.size  # 23.7 bytes

# 可以看到使用稀疏矩阵存储数据可以极大的减少数据的内存占用(23.7倍),提升运行的速度。

# 02 标准预处理工作流程
# 该步骤主要是进行质控以及选择细胞进行进一步分析。
# Seurat中一些常用的质控(QC)指标:
# 
# nFeature_RNA:每个细胞中检测到的基因的数量。
# 低质量的细胞或空液滴通常含有很少的基因;
# 细胞双胞体或多胞体可能表现出异常高的基因计数。
# 
# nCount_RNA:细胞内检测到的分子总数。即检测到的基因的总数量(一个基因可能会因为有多个转录本而被检测到多次),
# 因此分子总数与基因数量密切相关。
# 
# percent.mt:对应线粒体基因组的序列的百分比。低质量/濒死的细胞通常表现出广泛的线粒体污染。
# 我们通常使用 PercentageFeatureSet()函数来计算线粒体基因相关的reads的百分比。代码如下:


head(pbmc@meta.data)
row.names(pbmc)[grep("^MT-", rownames(pbmc)) ] #rownames = row.names

## ============= >>二选一
# (1)"[[]]"运算符可以向对象pbmc的元数据meta.data中添加列。meta.data是一个存储QC统计数据的好地方。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")  
# (2) 这是没有完全整理好的，没有MT开头，而是线粒体基因的通用名，总共只有13个基因
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "ATP6|ATP8|COX1|COX2|COX3|CYTB|ND1|ND2|ND3|ND4|ND4L|ND5|ND6")

VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA",  "percent.mt"),  ncol = 3)

# featureScatter是一种典型的可视化特征与特征之间关系的方法
plot1 <- FeatureScatter (pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter (pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

# 根据小提琴图和散点图的情况，个人认为可以选择如下的质控过滤指标，
# 即： nFeature_RNA > 200（最开始的min.features = 200）& nFeature_RNA < 2000 nCount_RNA < 7000 percent.mt < 7 
# 但是呢，为了和官方教程相统一，我们还是采用官方的过滤阈值，实际上相差并不多!
pbmc <- subset (pbmc, 
                subset = nFeature_RNA > 200  & 
                         nFeature_RNA < 2500 & 
                         percent.mt   < 5 )
pbmc 
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)
# 1 layer present: counts



# check data after QC
VlnPlot(pbmc, features =c("nFeature_RNA", "nCount_RNA",  "percent.mt"),  ncol = 3)
plot1 <- FeatureScatter (pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter (pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2


## 数据标准化 通过质控过滤完不需要的细胞之后呢，就要对数据进行标准化处理了！ LogNormalize：全局缩放的标准化方法（其实还有另一种标准化方法SCTransform，随后咱单独整理一篇！） 
# 为什么要做normalization呢？ 因为细胞间的异质性除了可能来自于细胞自身基因的差异表达，还有可能是由测序技术误差造成的。数据标准化的目的便是排除技术误差，让测序深度和文库大小不同的细胞的基因表达量具有可比性，

# 计算方法 we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), 
# and log-transforms the result. 
# log1p指的是log(x + 1)，
# 
pbmc <- NormalizeData (pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# Depreciated !! 结果储存在pbmc[["RNA"]]@data中,可以通过pbmc[["RNA"]]@data查看
# 使用 GetAssayData() 访问归一化后的数据 (通常是 log-normalized 数据)
head(GetAssayData(pbmc, assay = "RNA", slot = "data")[1:6, 1:25])
# 在 Seurat v5 (Assay5 对象) 中，内部数据的存储方式发生了变化，
# 主要的原始和归一化数据不再存储在 @data 和 @raw.data 槽位中，而是存储在 @layers 槽位中。
# 
# ✅ 解决方案：使用 GetAssayData() 函数
# 在任何版本的 Seurat 中，官方且推荐的访问 Assay 数据的安全方法是使用 GetAssayData() 函数，它可以兼容不同版本的内部存储结构。
# 如果您想访问原始计数 (Raw Counts)
head(GetAssayData(pbmc, assay = "RNA", slot = "counts")[1:6, 1:25])


# FindVariableFeatures鉴定高可变基因
pbmc <- FindVariableFeatures (pbmc, selection.method = "vst", nfeatures = 2000) 
# 鉴定10个最高变的基因
top10 <- head (VariableFeatures(pbmc), 10)
# 绘制高可变基因的火山图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot= plot1, points= top10, repel = TRUE)

combined_plot <- plot1 + plot2
combined_plot# << ========== plot is huge NOT visible in regular window
pdf("expression_and_variance.pdf", width = 12, height = 6)
#print(combined_plot) # << ========== plot is huge NOT visible in regular window
combined_plot# << ========== plot is huge NOT visible in regular window
dev.off()

# pdf(file=paste0(gene, ".PFS.pdf"), width=5.5, height=5, onefile=FALSE)
# print(surPlot)
# dev.off()

#　5 数据缩放（数据中心化） 
# 调整每个基因的表达量，使得细胞间的平均表达量为0 缩放每个基因的表达量，使得细胞间的方差为1 目的：这一步骤赋予每个基因在下游分析中同样的权重，
# 不至于使高表达基因占据主导地位 结果存储在pbmc[[“RNA”]]$scale.data中

all.genes <- rownames (pbmc)
# 以基因为单位对矩阵数据进行缩放
pbmc <- ScaleData (pbmc, features = all.genes) # features= scale.genes  will be faster  
pbmc <- RunPCA (pbmc, features = VariableFeatures (pbmc))

# PCA结果的可视化（其实不太重要）。 Seurat提供了几种方法来对PCA结果的细胞和基因进行可视化，包括VizDimReduce()、 DimPlot()和DimHeatmap()。 第一种：VizDimLoadings()：VizDimLoadings是一个R包中的函数，用于绘制维度加载图。维度加载图是一种用于可视化主成分分析（PCA）或其他降维方法中的变量与主成分之间的关系的图形。它显示了每个变量在每个主成分上的贡献程度，可以帮助我们理解变量在不同主成分上的权重和相关性。这种图形通常以散点图或箭头图的形式呈现。
VizDimLoadings (pbmc, dims = 1:2, reduction = "pca") #没啥用其实

DimPlot(pbmc, reduction = "pca") # #也没啥用,
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # 也没啥用
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # 也没啥用


ElbowPlot(pbmc) # 一般选择肘点附近（即本图中PC=10左右），作为数据集的PC值

#  fact from JackStraw
# JackStraw 分析是一种用于评估主成分分析中每个主成分显著性的方法，通过随机重采样的方式来模拟随机背景分布，然后将真实的主成分与随机背景进行比较，从而判断每个主成分是否包含真实的生物学信号。 pbmc：是一个 Seurat 对象，包含了单细胞 RNA 测序数据以及之前进行主成分分析（PCA）后的结果。 num.replicate = 100：指定随机重采样的次数，这里设置为 100 次。在每次重采样过程中，
# x = pvalue
pbmc=JackStraw(pbmc, num.replicate = 100)#: This runs the resampling procedure to generate the null distribution of gene scores.
pbmc=ScoreJackStraw(pbmc, dims = 1:15)
JackStrawPlot(pbmc,dim = 1:5)

#来自CSDN博主“生信小白要知道”的观点： 横坐标的PC是你接下来要选择分析的维度数，即多少组基因表达pattern。当选择一组的时候，即所有基因都在一个表达pattern里面（我们从生物学角度看也肯定是不可能的），我们很容易想到，所有细胞都在一个组内，肯定是存在很大差异的，所以方差值（方差=标准差（Standard Deviation）的平方）就会很大。当PC增多，即分组增多后，每个组内的差异也就慢慢变小了，方差值自然也就变小了。再往后，随着选的PC的增多，基本上区分不出来组内的差异了，这个时候方差基本就稳定了下来，此时在选取更大的PC值也就没有更大的意义了，并且还会增加计算压力。因此，个人认为关于PC选取的最合适原则是在拐点右侧3-5个维度即可（本数据拐点在7-8左右，我们选10刚刚好），宽松一点的话，可以设置在拐点右侧10以内。不过最重要的一点还是尽可能地多尝试，选取最适合自己的PC值，这一点也是官网教程所建议的。

#08 细胞聚类
# Seurat 应用了一种基于图表的聚类方法。重要的是,驱动聚类分析的距离度量(基于先前识别的PC) 保持不变。将细胞产制边界,然后将该生均中如K最近邻(KNN)图,在具有相似特征表达模式的细胞之间绘美的不同簇。
# 首先基于PCA空间中的欧几里得度量构建一个KNN图,然后基于局部邻域中的共享重叠(Jaccard 相似性)来精确任意两个细胞之间的边缘权重。这个步骤使用 FindNeighbors()函数执行,并将先前定义的前10个PC作为输入。

pbmc <- FindNeighbors(pbmc, dims = 1:10)
#为了将细胞聚类，我们接下来应用模块化优化技术，如Louvain算法(默认值)或SLM，迭代地将细胞聚类在一起。该步骤通过FindClusters()函数实现，其包含一个解析参数resolution，可以理解为就是分辨率
pbmc <- FindClusters(pbmc, resolution=0.5) # 0.4-1.2



# 假设您的 Seurat 对象名为 'pbmc'，并且您已经运行了 RunPCA(pbmc)

DimPlot(
    object = pbmc, 
    reduction = "pca",   # 指定使用 PCA 降维结果
    dims = c(1, 2)     ,  # 指定显示 PC1 (主成分 1) 和 PC2 (主成分 2)
    # group.by = "seurat_clusters", # 可选：根据聚类结果或自定义分组着色
    label = TRUE                  # 可选：显示聚类标签
)

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# 关于UMAP和TSNE图的展示简单做一个描述：可以说UMAP展现的距离更真实，TSNE图展示的结果似乎更美观（TSNE图呈现出来的结果中，各个细胞团都很紧促，最终呈现出来细胞都抱团分布，出来的图个人感觉更美观一丢丢）。
pbmc<-RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

# 计算所有簇的markers Seurat教程仅设置了一个参数，就是only.pos，该参数的意思是只保留正标记（positive markers），也就是只保留相对于其他簇而言，目标簇显著高表达的基因标记。 FindAllMarkers() 函数：这是 Seurat 包中的一个函数，用于对 Seurat 对象 pbmc 中的所有细胞聚类进行差异表达分析，找出每个聚类的标记基因。 only.pos = TRUE：该参数指定只返回在特定聚类中高表达的基因，即正差异表达的基因，忽略那些在其他聚类中高表达的基因。 -min.pct = 0.25：是一个过滤条件，要求目标基因至少在 25% 的当前聚类细胞或至少 25% 的其他聚类细胞中表达，才会被纳入差异表达分析，这样可以减少在极少数细胞中表达的基因对分析结果的干扰。 logfc.threshold = 0.25：表示只考虑平均对数 2 倍变化值大于 0.25 的基因，即该基因在当前聚类中的表达量至少是其他聚类的 2^0.25 = 1.19 倍，这有助于筛选出表达差异更显著的基因。


# 找出每个细胞簇的标记物,与所有剩余的细胞进行比较,只报告阳性细胞
#pbmc.markers<-FindAllMarkers (pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.25)
pbmc.markers<-FindAllMarkers (pbmc, only.pos = TRUE)

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


#寻找cluster2的所有其他簇的 markers
cluster2.markers <- FindMarkers(pbmc, ident.1 =2, min.pct = 0.25)
#寻找cluster5中与cluster0和cluster3 不同的所有markers
cluster5.markers <- FindMarkers(pbmc, ident.1=5, ident.2 = c(0, 3) , min.pct = 0.25)
                            

# 实际上，到了这里单细胞Seurat流程的分析基本上走完了，在这里保存RDS是比较合适的，因为里面也存好了计算的差异基因的情况，到时候重新加载RDS文件时，不需要运算即可获取。同时与注释获得完整的细胞类型还有一定的时间要求。所以，我们选择在这里保存第一版分析数据的RDS。
#保存结果为rds文件

output_file <- "../output/pbmc_tutorial.rds"

# 2. 从完整路径中提取文件夹路径
output_dir <- dirname(output_file)

# 3. 检查文件夹是否存在，如果不存在则创建

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("已创建文件夹:", output_dir, "\n")
} else {
    cat("文件夹已存在:", output_dir, "\n")
}

# 4. 保存 Seurat 对象
saveRDS(pbmc, file = output_file)
cat("Seurat 对象已成功保存到:", output_file, "\n")


#11 基因表达可视化 这一部分是专门添加进来的内容，特意和第10步的内容相区别。关于基因表达可视化，俺有五种方法，任君挑选。至于图的含义，我就不做阐述了，已经挺直观啦！And，后续我们也会出基因表达可视化绘图的内容，可以期待一下
# 第1种：VizDimLoadings()
# 第2种：DimPlot()：（PCA, t-SNE, UMAP）
# 第3种：DimHeatmap()
# 第4种：VlnPlot()
# #Dicking joint handwidth of A 319
# 第5种：山脊图
     
# 按簇展示特定基因表达水平
VlnPlot (pbmc, features = c("MS4A1", "CD79A")) 


# X 取离散整数值时, log(X+1) 也会取固定的、离散的值：
# UMI 原始计数 XLog 转换值 log(X+1) (以 e 为底)为什么出现直线
# log(1)  大量的细胞聚集在 Y=0 的线上 
# log(2) 所有 UMI=1 的细胞聚集在 Y≈0.693 的线上
# log(3)1.099$所有 UMI=2 的细胞聚集在 Y≈1.099 的线上
# log(4)1.386$所有 UMI=3 的细胞聚集在 Y≈1.386 的线上
# log(11)2.398$所有 UMI=10 的细胞聚集在 Y≈2.398 的线上
# log(101)4.615$所有 UMI=100 的细胞聚集在 Y≈4.615 的线上
# 
# #可以绘制行数据
VlnPlot (pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


features <- c("LYZ", "CCL5", "IL32", "PTPRCAP")
RidgePlot(pbmc, features = features, ncol = 2)
#Dicking joint handwidth of A 319

FeaturePlot (pbmc, features = features) # 我命名其为基因表达点图

# 点阵图
DotPlot (pbmc, features = features) + RotatedAxis()

DoHeatmap(subset(pbmc, downsample =100), features = features,  size = 3)
#DoHeatmap()为给定的细胞和特征生成一个表达热图。在这种情况下,我们为每个集群绘制前20个标记(或所有标记,如果小于20)。
pbmc.markers%>%
    group_by (cluster)%>%
    top_n( n = 100, wt=avg_log2FC) -> top100

top100
DoHeatmap (pbmc, features = top100$gene) + NoLegend()

#这里记住一点，基因表达量的图是用来显示基因在各簇中的表达情况的，这个为该基因是否可以作为本簇的marker提供了直接的证据来源，如和已有报道相互印证（如XXX细胞特异性表达该基因），则可通过该基因对某一簇细胞进行细胞类型注释。细胞类型注释是一个人工主观的过程（自动注释的应用范围很少，只聚焦在极少数物种中，大多数的细胞类型注释还是依赖人工进行的），这一过程已游离于上述前10步分析之外，所以这就是为什么我在第10步结束的时候保存.rds的原因，不过这里也确实没有细讲如何进行细胞类型注释，咱们后续专门出一期和大家唠唠细胞注释的内容！ 12 细胞簇命名 假定我们走完上述11步，已经完成了细胞类型注释工作，即知道了每个簇对应的细胞类型。那么现在我们要如何将这些簇更改成对应的细胞名称呢？看我手法：
# 如果您的聚类 ID 数量很多，建议使用 data.frame 或 Excel 表格来整理映射关系，然后导入 R 中转换为命名向量。
#把每个统对应的新名字写成一个能合
#这里注意,如果有两个候注释到的名字一样,
#要在对应数字位置都写一遍<<==即一定保证这个集合的元素数和保的个数是一致的8
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Мемогу CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "platelet")

new.cluster.ids

#这一步是把pbmc的levels传输给new.cluster.ids作为他的name属性
# 1. 创建命名向量
# 格式为：New Name = Old ID 
# 注意：在 R 中创建命名向量时，通常是将 '值' 赋给 '名称'
new.cluster.ids <- c("CD4 T cells", "B cells", "CD8 T cells", 
                     "NK cells", "Monocytes") 
new.cluster.ids
names(new.cluster.ids) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8") 

# 2. 执行重命名操作
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 3. 验证结果
head(Idents(pbmc)) # 查看前几个细胞的新 ID


3 #最后,再展示一下注释后的细胞图谱
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size= 0.5) + NoLegend()














# 到这里，Seurat的基础分析算全部完成了，这时候我们可以保存一下最完整的信息。
saveRDS (pbmc, file="../output/pbmc3k_final.rds")