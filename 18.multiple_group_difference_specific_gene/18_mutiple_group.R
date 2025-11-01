setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/18.multiple_group_difference")
library(limma)
library(ggpubr)

rt1=read.table("A.txt",sep="\t",header=T,check.names=F)
rt2=read.table("B.txt",sep="\t",header=T,check.names=F)
rt3=read.table("C.txt",sep="\t",header=T,check.names=F)
rt4=read.table("D.txt",sep="\t",header=T,check.names=F)
# rt1=as.matrix(rt1)
# rt2=as.matrix(rt2)
# rt3=as.matrix(rt3)
# rt4=as.matrix(rt4)
rownames(rt1)=rt1[,1]
rownames(rt2)=rt2[,1]
rownames(rt3)=rt3[,1]
rownames(rt4)=rt4[,1]
rt1=rt1[,2:ncol(rt1)]
rt2=rt2[,2:ncol(rt2)]
rt3=rt3[,2:ncol(rt3)]
rt4=rt4[,2:ncol(rt4)]

bindrt=cbind(rt1,rt2,rt3,rt4)
dimnames=list(rownames(bindrt),colnames(bindrt))
data=matrix(as.numeric(as.matrix(bindrt)),nrow=nrow(bindrt),dimnames=dimnames)
write.table(data, file="all.txt", sep="\t", quote=F, row.names=T)

ps=sapply(strsplit(colnames(data),""), "[", 1)
Type=as.matrix(groups)
colnames(Type)="group"
ids=as.matrix(colnames(as.matrix(data)))
colnames(ids)="id"
groupbind=cbind(ids, groups)
groupbind=as.matrix(groupbind)
rownames(groupbind)=groupbind[,1]
groupbind=as.data.frame(groupbind[,2:ncol(groupbind)])
colnames(groupbind)="group"
write.table(groupbind, file="groups.txt", sep="\t", quote=F, row.names=T)

####ample=intersect(colnames(data), row.names(groupbind))
data=t(data[drop=F,,samSample])
cli=groupbind[samSample,,drop=F]
rt=cbind(data, cli)

####="TSPAN6"
Types="group"
data=rt[c(gene, Types)]
colnames(data)=c(gene, "Types")
data=data[(data[,"Types"]!="unknow"),]
#???up=levels(factor(data$Types))
data$Types=factor(data$Types, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#???plot=ggboxplot(data, x="Types", y=gene, fill="Types",
                  xlab=Types,
                  ylab=paste(gene, " expression"),
                  legend.title=Types)+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
#???(file=paste0("TypesCor_", Types, ".pdf"), width=5.5, height=5)
print(boxplot)
dev.off()
