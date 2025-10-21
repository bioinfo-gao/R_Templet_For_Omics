#biocLite("limma")
#install.packages("gplots")
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

###limma包分析TCGA数据要用FPKM，不然矫正后的P值会很大
###limma包分析TCGA数据要用FPKM，不然矫正后的P值会很大
###limma包分析TCGA数据要用FPKM，不然矫正后的P值会很大
###limma包分析TCGA数据要用FPKM，不然矫正后的P值会很大
#setwd("")
rt=read.table("combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)
#多行取平均值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data2=as.data.frame(data)

#以01A和11A分组，正常放前面，肿瘤放后面
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp_data_N) 
rt= cbind(exp_data_N, exp_data_T)
#write.table(rt,file = "groupout.txt",sep="\t",quote=F)



#加载limma包，用于校正和比较差异
#rt=read.table("groupout.txt",sep="\t",header=T,check.names=F,row.names = 1)

#normalize以下为对数据的校正
pdf(file="rawBox.pdf")
#画箱线图
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
#校正
rt=normalizeBetweenArrays(rt)
pdf(file="normalBox.pdf")
#画箱线图
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

group1=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       #正常组样品数目
treatNum=length(group1[group1==0])     #肿瘤组样品数目
#differential差异分析
class <- c(rep("con",conNum),rep("treat",treatNum))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")
#算方差
df.fit <- lmFit(rt,design)
df.matrix<- makeContrasts(con - treat,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',n=Inf) 
#写入表格
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#找出差异两倍以上，pvalue小于0.05，写入表格
diffLab <- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)

#差异基因表达水平，用于共表达
diffExpLevel <- rt[rownames(diffLab),]
write.table(diffExpLevel,file="diffExpLevel.xls",sep="\t",quote=F)

#可视化
gene="THBS2" 
data=t(rt[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp$gene=log2(exp$gene+1)
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                  xlab="",
                  ylab=paste0(gene, " expression"),
                  legend.title="Type",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
#输出图片
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

