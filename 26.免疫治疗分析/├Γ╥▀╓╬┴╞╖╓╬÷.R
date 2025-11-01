
#引用包
library(limma)
library(ggpubr)

tciaFile="TCIA.txt"        #免疫打分文件
expFile="geneExp.txt"      #表达数据文件

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#删掉正常样品
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#根据目标基因表达量对样品进行分组
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)

#读取TCIA的打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(ips), row.names(data))
ips=ips[sameSample, , drop=F]
data=data[sameSample, "Type", drop=F]
data=cbind(ips, data)

#设置比较组
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c("Low", "High"))
group=levels(factor(data$Type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#可视化
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "Type")]
	gg1=ggviolin(rt, x="Type", y=i, fill = "Type", 
	         xlab="", ylab=i,
	         legend.title=gene,
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=4.8, height=4.25)
	print(gg1)
	dev.off()
}



