# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")                 # #biocLite("limma")
# install.packages(c( "gplots", "beepr", "tidyverse" ))
#install.packages("psych")  # for data analysis # install.packages("psycho") # focus on psychcology
# install.packages("here")
#library(here) # here 包来构建跨平台的路径，它会自动处理分隔符问题。

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
library(psych)

getwd()
setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\02.差异分析//")
    
# 2 min to read 
rt=read.table("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/01.New_TCGA/combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)

dim(rt)
class(rt)
head(rt)
rt[ 1:3, 1:5 ]

rt=as.matrix(rt)

rownames(rt) # rownames(rt)=rt[,1]

exp=rt[ ,3:ncol(rt)] # exp=rt[ ,2:ncol(rt)]

exp[1:2, 1:5]
class(exp) # [1] "data.frame"

dimnames=list(rownames(exp), colnames(exp))

dimnames # gene names & sample names         ?nrow

# temp1 =  (as.matrix(exp) )
# head(temp1)
# class(temp1)
# describe(temp1)
# psych::describe(exp) # psych::describe(temp1)
# summary(temp1) 
# ?summary
# temp1[1:2, 1:6]

# data=matrix( as.numeric( as.matrix(exp) ),  nrow=nrow(exp),   dimnames = dimnames) # nrow=nrow(exp) <<== is NEEDED !!
# data=matrix( as.numeric( as.matrix(exp) ),  nrow=nrow(exp),   dimnames = dimnames) # nrow=nrow(exp) <<== is NEEDED !!
# 
# data[1:3, 1:3]

# Condense a microarray data object so that values for within-array replicate probes 
# are replaced with their average.
#?avereps()
# Not needed at all here for NGS
# data=avereps(data) 
data = exp
data=data[rowMeans(data)>0,]
colnames(data)

# data2=as.data.frame(data)
# 
# data[1:3, 1:3]
# data2[1:3, 1:3]
# colnames(data2)

#??01A??11A???飬??????ǰ?棬?????ź???

exp_data_T = data %>% dplyr::select(str_which( colnames(.),  "-01A")) #%>% can NOT be replaced !! 

exp_data_T[1:5, 1:5] 
nT = ncol(exp_data_T) 

nT

exp_data_N = data %>% dplyr::select(str_which(colnames(.), "-11A")) # %>% can NOT be replaced 

exp_data_N[1:5, 1:5] 
nN = ncol(exp_data_N) 

nN

rt= cbind(exp_data_N, exp_data_T)


dim(rt)

rt[1:2,1:5]

getwd()

pdf(file="rawBox.pdf", width = 100, height = 5) #inches, always take care of "," in Chinese
boxplot(rt, col = "blue", xaxt = "n", outline = F) 
dev.off()


# normalizede reads  
rt=normalizeBetweenArrays(rt)

pdf(file="normalBox.pdf", width = 100, height = 5)
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

# 01 或 02 等代表肿瘤（Tumor）样本。
# 10 或 11 等代表正常（Normal）样本。
group1 = sapply(strsplit(colnames(rt),"\\-"), "[", 4) # 提取每个子向量的第四个元素。
#  使用空字符串 ("") 作为分隔符，将 group1 中的每个字符串分割成单个字符。
group1 = sapply(strsplit(group1      ,""   ), "[", 1) # 这里的意思是逐个字符分隔string ，所以不需要分隔符
group1 = gsub("2", "1", group1) # 在 将所有出现的字符串 "2" 替换为字符串 "1"。
group1

conNum  =length(group1[group1==1])     
treatNum=length(group1[group1==0])     

#differential
class <- c(rep("con",conNum),rep("treat",treatNum))  
class

# factor() 确保R将这个变量视为一个分类变量（组别）
# + 0: 这是最关键的部分，它表示不包含截距项（Intercept）。
# 在 limma 分析中，强烈建议使用 + 0 的细胞平均数模型
#（即设计矩阵中包含所有组的列）。这样做的好处是：
#你可以直接定义任何两组之间的差异（例如 B - A），而不需要考虑基线参照组，
design <- model.matrix( ~factor(class) +0 ) 

colnames(design)  <-  c("con","treat")
colnames(design)
design


# ==============================================================
df.fit    <- lmFit(rt, design)
df.matrix <- makeContrasts(con - treat,levels=design)
fit       <- contrasts.fit(df.fit,df.matrix)

?eBayes # compute moderated t-statistics, moderated F-statistic, and log-odds of differential
# expression by empirical Bayes moderation of the standard errors towards a global value.
fit2 <- eBayes(fit)


allDiff = topTable(fit2, adjust='fdr', n=Inf )  # n=Inf	(数量)	要返回的基因数量


write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)


# abs fold > 2 AND FDR < 0.05
diff<- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
writdiffLab

e.table(diffLab,file="diffExp.xls",sep="\t",quote=F)

ExpLevel <- rt[rownames(diffLab),]
write.table(diffExpLevel,file="diffExpLevel.xls",sep="\t",quote=F)


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
    stat_compare_means(comparisons=my_comparisons,
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                       label = "p.signif")
nt(boxplot)
dev.off()

