

#BiocManager::install("edgeR")
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
library(future.apply)
foldChange=2
padj=0.05


rt=read.table("combined_RNAseq_counts.txt",sep="\t",header=T,check.names=F) #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)


#以01A和11A分组，正常放前面，肿瘤放后面
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 

exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 

data= cbind(exp_data_N, exp_data_T)


group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)

conNum  =length(group1[group1==1])       #正常组样品数目
treatNum=length(group1[group1==0])     #肿瘤组样品数目


group=c(rep("normal",conNum),rep("tumor",treatNum)) #按照自己的数据更改正常组和肿瘤组的数量
design <- model.matrix(~group)

# ??DGEList

y <- DGEList(counts=data,group=group)#构建列表 edgeR

keep.exprs<-filterByExpr(y,group = group)

y <- y[keep.exprs,,keep.lib.sizes=FALSE]

dim(y)

nsamples<-ncol(y)
col<-brewer.pal(nsamples,"Paired")
y <- calcNormFactors(y)#计算样本内标准化因子
y[["samples"]][["norm.factors"]]

y2 <- y
y2$samples$norm.factors <- 1
y2$counts[,1] <- ceiling(y2$counts[,1]*0.05)
y2$counts[,2] <- y2$counts[,2]*5

#标准化前后箱线图
par(mfrow=c(1,2))
lcpm1<-cpm(y2,log=TRUE)

pdf(file="rawBox.pdf",width=20,height=20)
boxplot(lcpm1,las=2,col=col,main="",xaxt = "n")
title(main="A.Example:Unnormalised data",ylab="Log-cpm")
dev.off()

y2<-calcNormFactors(y2)
y2$samples$norm.factors

lcpm2<-cpm(y2,log=TRUE)

pdf(file="normalBox.pdf",width=20,height=20)
boxplot(lcpm2,las=2,col=col,main="",xaxt = "n")
title(main="B.Example:Normalised data",ylab="Log-cpm")
dev.off()

y <- estimateCommonDisp(y2)#计算普通的离散度
y <- estimateTagwiseDisp(y)#计算基因或miRNA范围内的离散度

et <- exactTest(y,pair = c("normal","tumor"))#进行精确检验

topTags(et)#输出排名靠前的差异miRNA信息

ordered_tags <- topTags(et, n=100000)#将差异信息存入列表

#剔除FDR值为NA的行

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff

newData=y$pseudo.counts #将y中的counts信息通过管道存入newData

#将差异miRNA的logF、logCPM、PValue、FDR存入表格
write.table(diff,file="edgerOut.xls",sep="\t",quote=F)

#将差异miRNA：p值大于规定值、logFC大于规定值的logF、logCPM、PValue、FDR存入表格
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)

#将上调miRNA的logF、logCPM、PValue、FDR存入表格
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)

#将下调miRNA的logF、logCPM、PValue、FDR存入表格
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

#将差异miRNA在tumor和normal样本中的表达值存入表格
#
normalizeExp=rbind(id=colnames(newData),newData)

write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F) #输出所有基因校正后的表达值（normalizeExp.txt）

diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)#输出差异基因校正后的表达值（diffmRNAExp.txt）


#可视化
#样品分组
gene="DCAF5"
data=t(newData[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
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
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)

print(boxplot)

dev.off()



