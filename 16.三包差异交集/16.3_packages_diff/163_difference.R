library(edgeR)
library(data.table)
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)
library("DESeq2")
library(VennDiagram)
padj = 0.05
foldChange= 2
#####limma分析#####
rt1=read.table("combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)
#多行取平均值
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
data1=avereps(data1)
data1=data1[rowMeans(data1)>0,]
data1=as.data.frame(data1)
#以01A和11A分组，正常放前面，肿瘤放后面
exp1_data_T = data1%>% dplyr::select(str_which(colnames(.), "-01A")) # 匹配列名或用下示写法
nT = ncol(exp1_data_T) 
exp1_data_N = data1%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp1_data_N) 
rt1= cbind(exp1_data_N, exp1_data_T)
#校正
rt1=normalizeBetweenArrays(rt1)
group1=sapply(strsplit(colnames(rt1),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum1=length(group1[group1==1])       #正常组样品数目
treatNum1=length(group1[group1==0])     #肿瘤组样品数目
#differential差异分析
class <- c(rep("con",conNum1),rep("treat",treatNum1))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")
#算方差
df.fit <- lmFit(rt1,design)
df.matrix<- makeContrasts(con - treat,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDEG1 = topTable(fit2,coef=1,n=Inf,adjust="BH") 
allDEG1 = na.omit(allDEG1)
#提取基因差异显著的差异矩阵
padj = 0.05
foldChange= 2
diff_signif1 = allDEG1[(allDEG1$adj.P.Val < padj & 
                          (allDEG1$logFC>foldChange | allDEG1$logFC<(-foldChange))),]
diff_signif1 = diff_signif1[order(diff_signif1$logFC),]
save(diff_signif1, file = 'limma_diff1.Rdata')



#####DESeq2差异分析#####
rt <- read.table( "combined_RNAseq_counts.txt",header=T,sep="\t",comment.char="",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)


exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)


group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       #正常组样品数
treatNum=length(group1[group1==0])     #肿瘤组样品数

count <- floor(data)#下取
#count <- ceiling(count)#上取
# 预处理，过滤低丰度的数据
countData <- count[apply(count, 1, sum) > 0 , ]
# 读取样本分组信息
data=colnames(countData)
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("id", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp=as.matrix(exp)
rownames(exp)=exp[,1]
exp=exp[,2:ncol(exp)]
exp=as.data.frame(exp)
colnames(exp)=c("condition")
write.table(exp, file="group.txt",sep="\t",quote=F)
colData <- read.table( "group.txt",header=T,sep="\t",row.names=1,comment.char="",check.names=F)
# 构建DESeq2中的对象
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
# 指定哪一组作为对照组
dds$condition <- relevel(dds$condition, ref = "Normal")
#计算每个样本的归一化系数
dds <- estimateSizeFactors(dds)
#估计基因的离散度
dds <- estimateDispersions(dds)
#差异分析
dds <- nbinomWaldTest(dds)
dds <- DESeq(dds) 
allDEG2 <- as.data.frame(results(dds))
#提取基因差异显著的差异矩阵
padj = 0.05
foldChange= 2
diff_signif2 = allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif2 = diff_signif2[order(diff_signif2$log2FoldChange),]
save(diff_signif2, file = 'DESeq2_diff2.Rdata')


#####degeR分析#####
rt=read.table("combined_RNAseq_counts.txt",sep="\t",header=T,check.names=F) #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)
group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       #正常组样品数
treatNum=length(group1[group1==0])     #肿瘤组样品数
group=c(rep("normal",conNum),rep("tumor",treatNum)) #按照自己的数据更改正常组和肿瘤组的数???
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)#构建列表
y <- calcNormFactors(y)#计算样本内标准化因子
y <- estimateCommonDisp(y)#计算普通的离散
y <- estimateTagwiseDisp(y)#计算基因或miRNA范围内的离散
et <- exactTest(y,pair = c("normal","tumor"))#进行精确检
topTags(et)#输出排名靠前的差异miRNA信息
ordered_tags <- topTags(et, n=100000)#将差异信息存入列
#剔除FDR值为NA的行
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
#提取基因差异显著的差异矩阵
padj = 0.05
foldChange= 1
diff_signif = allDiff[(allDiff$FDR < padj & (allDiff$logFC>foldChange | allDiff$logFC<(-foldChange))),]
diff_signif = diff_signif[order(diff_signif$logFC),]
save(diff_signif, file = 'edger_diff.Rdata')


#####可视化#####
edgeR = rownames(diff_signif)
dim(diff_signif)
limma = rownames(diff_signif1)
dim(diff_signif1)
DESeq2 = rownames(diff_signif2)
dim(diff_signif2)

venn.diagram(
  x = list(
    'edgeR' = edgeR,
    'limma' = limma,
    'DESeq2' = DESeq2
  ),
  filename = 'VN.png',
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main = "三种包的差异表达基因比较",
  main.cex = 1.2
)

#####保存表格#####
edgeR=as.data.frame(edgeR)
limma=as.data.frame(limma)
DESeq2=as.data.frame(DESeq2)
sameSample=intersect(edgeR$edgeR, limma$limma)
sameSample=as.data.frame(sameSample)
sameSample=intersect(DESeq2$DESeq2,sameSample$sameSample)
#保存交集基因
sameSample1=as.data.frame(sameSample)
write.table(sameSample1, file="merge_genes.xls",sep="\t",quote=F)
#提取counts
data1=exp[sameSample,,drop=F]
write.table(data1, file="merge_counts.xls",sep="\t",quote=F)
#提取FPKM
data2=exp1[sameSample,,drop=F]
write.table(data2, file="merge_FPKM.xls",sep="\t",quote=F)
