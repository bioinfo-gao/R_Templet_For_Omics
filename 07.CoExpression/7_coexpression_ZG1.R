#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")
install.packages("future")
library(future)


library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

gene="ARL15"               # for plot 
corFilter=0.3           
pFilter=0.05            

# the default \ is not allowed in R, must use Either \\ or / == ZG
data_dir = "C:/Users/zhen-/Code/R_code/R_For_DS_Omics/01.New_TCGA"
#data_dir = "C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\01.New_TCGA"

expFile= file.path(data_dir,  "combined_RNAseq_FPKM.txt"  )    
expFile

rt=rt1=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rt[1:2, 1:5]

rownames(rt)=make.unique(rt[,1])  # ZG

exp=rt[,3:ncol(rt)] # not 2
exp[1:2, 1:5]

dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data[1:2, 1:6]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=log2(data+1)
data[1:2, 1:6]



#??ȡĿ??????????��
x=as.numeric(data[gene,])
#?Ի???????ѭ?????????????Լ???
res3 <- data.frame()
genes <- colnames(rt2)[-c(1:2)]

plan(multisession) # library(future) 
# launch a number of background R sessions equal to the number of available cores detected by availableCores()

system.time( res3 <- future_lapply(rownames(data)), function(j){ # future lib 
    print(j)
    if(gene==j){next}
    y=as.numeric(data[j,])
    corT=cor.test(x, y, method = 'pearson')
    cor=corT$estimate
    pvalue=corT$p.value
    outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
    
    if((abs(cor)>corFilter) & (pvalue<pFilter)){
        
        df1=as.data.frame(cbind(x,y))
        p1=ggplot(df1, aes(x, y)) + 
            xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
            geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
            stat_cor(method = 'pearson', aes(x =x, y =y))
        pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
        print(p1)
        dev.off()   }
    
} 
)



res3 <- data.frame(do.call(rbind,res3))
for(j in rownames(data)){
    if(gene==j){next}
    y=as.numeric(data[j,])
    corT=cor.test(x, y, method = 'pearson')
    cor=corT$estimate
    pvalue=corT$p.value
    outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
    #?????????????Ļ???
    if((abs(cor)>corFilter) & (pvalue<pFilter)){
        #???ӻ?
        df1=as.data.frame(cbind(x,y))
        p1=ggplot(df1, aes(x, y)) + 
            xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
            geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
            stat_cor(method = 'pearson', aes(x =x, y =y))
        pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
        print(p1)
        dev.off()
    
#?????????Խ????ļ?
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
outTab=outTab[abs(as.numeric(outTab$cor))>corFilter & as.numeric(outTab$pvalue)<pFilter,]
write.table(file="corSig.txt", outTab, sep="\t", quote=F, row.names=F)



