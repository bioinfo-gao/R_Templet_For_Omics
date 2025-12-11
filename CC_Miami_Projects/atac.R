##################################################################################################
#     Profiling accessible chromatin landscapes in mouse tissues by ATAC-seq R code              #
##################################################################################################

library(DESeq2)
library(viridis)
library(reshape2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(heatmaply)
library(webshot)
library(Rtsne)
library(RColorBrewer)

filecount = "chromatin.accessibility.raw.count.txt"
fileTF = "TF.motif.enrichment.txt"
outpath = "Users/out/"

##-------------------------------------------------------------##
##              Normalize raw readscount by RPM                ##
##-------------------------------------------------------------##

data=read.table(filecount,header = T)
row.names(data)=data[,1]
data=data[,-1]

colsum=colSums(data)
RPM<-function(x){return(x/(colsum/1000000))}
data_RPM=as.data.frame(t(apply(data,1,RPM)))

##-------------------------------------------------------------##
##                     Correlation test                        ##
##-------------------------------------------------------------##

##The correlation between all samples

col1 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","white","yellow","#FF7F00","red","#7F0000"))
heatmaply(cor(log10(data_RPM[,-c(1:12)]+1)),colors = col1(100),column_text_angle = 45,file="cor_heatmap.html")
webshot::webshot("cor_heatmap.html", file="cor_heatmap.pdf", delay=2)

##The correlation between technical replicates or biological replicates

Lab.palette <- colorRampPalette(c("white","sky blue","blue","royal blue","dodger blue","medium sea green","orange", "yellow"),space="rgb")

#Cerebrum
smoothScatter(log10(data_RPM[,c(21,22)]+1),colramp = Lab.palette,
              main= paste("cor = ",cor(log10(data_RPM[,21]+1),log10(data_RPM[,22]+1)),sep = ""))
smoothScatter(log10(data_RPM[,c(21,53)]+1),colramp = Lab.palette,
              main= paste("cor = ",cor(log10(data_RPM[,21]+1),log10(data_RPM[,53]+1)),sep = ""))

##The correlation between our data and published studies

#Heart
smoothScatter(log10(data_RPM[,c(3,23)]+1),colramp = Lab.palette,
              main= paste("cor = ",cor(log10(data_RPM[,3]+1),log10(data_RPM[,23]+1)),sep = ""))
#Lung
smoothScatter(log10(data_RPM[,c(10,29)]+1),colramp = Lab.palette,
              main= paste("cor = ",cor(log10(data_RPM[,10]+1),log10(data_RPM[,29]+1)),sep = ""))

##-------------------------------------------------------------##
##               Identify tissue specific peaks                ##
##-------------------------------------------------------------##

##calculate average between replicates

nor_RPM=data_RPM[,-c(1:12)]

test=as.data.frame(melt(nor_RPM[1,]))
test$group=sapply(strsplit(as.character(test$variable),'_'),"[", 2)
name=sapply(strsplit(as.character(test$variable),'_'),"[", 2)
mean=as.data.frame(t(tapply(test$value,test$group,mean)))

out=apply(nor_RPM,1,function(x){
    test=as.data.frame(melt(x))
    test$group=name
    mean=t(tapply(test$value,test$group,mean))
    return(mean)
})

ave=as.data.frame(t(out))
colnames(ave)=colnames(mean)

##Using entropy method to identify tissue specific peaks

ratio=as.data.frame(t(apply(ave,1,function(x) x/sum(x))))
entropy <- function(a){
    b<- -sum(a*log2(a))
    return(b)
}
entro=as.data.frame(apply(ratio,1,entropy))

##distribution of entropy scores

omit_entro=as.data.frame(entro[-which(entro[,1]=="NaN"),])
colnames(omit_entro)=c("entropy")
sta=as.data.frame(seq(from=1.3,to=4.4,by=0.01))

for(i in 1:nrow(sta)){
    
    sub_entropy=subset(omit_entro,omit_entro$entropy<sta[i,1])
    sta[i,2]=nrow(sub_entropy)
    
}

colnames(sta)=c("entropy","peak_num")

ggplot(data = sta, mapping = aes(x = entropy, y = peak_num)) + 
    geom_point(size=0.3)+
    theme_bw()+
    xlab("entropy score")+
    ylab("peak number")

##plot heatmap

plot=cbind(ave,entro)

sub_plot=subset(plot,plot$`apply(ratio, 1, entropy)`<3.5)#Peaks with entropy less than 3.5 were selected as tissue-restricted peaks.
order_sub_plot=sub_plot[order(sub_plot$`apply(ratio, 1, entropy)`),]
tissue_name=colnames(sub_plot[,-21])
for(i in seq(nrow(order_sub_plot))){
    a=as.matrix(order_sub_plot[i,-21])
    order_sub_plot[i,22]=tissue_name[which(a==a[which.max(a)],arr.ind=T)[2]]
}
re_order_sub_plot=order_sub_plot[order(order_sub_plot$V22),]
re_order_sub_plot=re_order_sub_plot[,c(1,2,3:11,13,14,12,15,16,18,19,17,20)]#sort tissues manually
pheatmap(log2(re_order_sub_plot+1),cluster_rows = FALSE,cluster_cols = FALSE,color= inferno(50),main="Tissue specific peaks",show_rownames = F)

##subset tissue specific peaks
forsub=cbind(ratio,entro)
sub=subset(forsub,forsub$`apply(ratio, 1, entropy)`<3.5)
tissue_name=colnames(sub[,-21])
for(i in seq(nrow(order_sub_plot))){
    a=as.matrix(sub[i,-21])
    b=i
    for(i in 1:20){
        if(sort(a,decreasing = T)[i] > 0.1){sub[b,21+i] = tissue_name[which(a==sort(a,decreasing = T)[i])] }
        else{sub[b,21+i] = "FALSE"}
    }
}
nth_1=subset(sub[22],sub[22] != "FALSE")#subset peaks in one tissue
nth_2=subset(sub[23],sub[23] != "FALSE")#subset peaks in two tissues
nth_3=subset(sub[24],sub[24] != "FALSE")#subset peaks in three tissues
nth_4=subset(sub[25],sub[25] != "FALSE")#subset peaks in four tissues
nth_5=subset(sub[26],sub[26] != "FALSE")#subset peaks in five tissues
nth_6=subset(sub[27],sub[27] != "FALSE")#subset peaks in six tissues
colnames(nth_1)=c("tissue")
colnames(nth_2)=c("tissue")
colnames(nth_3)=c("tissue")
colnames(nth_4)=c("tissue")
colnames(nth_5)=c("tissue")
colnames(nth_6)=c("tissue")
nth_1$id=row.names(nth_1)
nth_2$id=row.names(nth_2)
nth_3$id=row.names(nth_3)
nth_4$id=row.names(nth_4)
nth_5$id=row.names(nth_5)
nth_6$id=row.names(nth_6)
ALL_peaks=rbind(nth_1,nth_2,nth_3,nth_4,nth_5,nth_6)

for(i in 1:20){
    test=subset(ALL_peaks,ALL_peaks$tissue==tissue_name[i])[2]
    write.table(test,paste0(outpath,"entropy_3.5_",tissue_name[i],"_peaks.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}

##-------------------------------------------------------------##
##                           t-SNE analysis                    ##
##-------------------------------------------------------------##

sub_data_RPM=data_RPM[,-c(1:12)]
sub_data_RPM$ID=rownames(sub_data_RPM)
ALL_peaks$ID=row.names(ALL_peaks)

use=merge(sub_data_RPM,ALL_peaks,by="ID")
use=use[,-c(1,68,69)]
tissue=as.data.frame(colnames(use))

a=sub("Female_","",tissue[,1])
b=sub("Male_","",a)
c=sub("_Rep1","",b)
d=sub("_Rep2","",c)

train<- as.data.frame(t(use))

train$label=d
train=train[,c(26197,1:26196)]
train$label<-as.factor(train$label)

tsne<- Rtsne(as.matrix(train[,-1]), dims = 2, perplexity=5, verbose=TRUE, max_iter = 1500)

coor=as.data.frame(tsne$Y)
coor$Tissue=d

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colnames(coor)=c("tSNE1","tSNE2","Tissue")
shape=as.integer(c(0,18,5,16,15,14,13,12,6,10,9,8,7,6,5,4,3,2,1,0))
ggplot(data=coor,aes(x=tSNE1, y = tSNE2,colour=Tissue,shape=Tissue))+
    geom_point(size=3) +
    theme_gray() +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
    scale_shape_manual(values = shape)+
    scale_color_manual(values = getPalette(20))+
    ggtitle("ATAC-Tissue-Correlation")

##-------------------------------------------------------------##
##               Identify tissue specific TFs                  ##
##-------------------------------------------------------------##
#calculate -log10 P-valus and CV
data=read.table(fileTF,header = T,sep = "\t")
data[data==0]<- 1
row.names(data)=data[,1]
data=data[,-1]
log_data=-log10(data[,-1])
log_data$Motif=data$family.motif
log_data$sd=apply(log_data[,-21],1,sd)
log_data$mean=apply(log_data[,-21],1,mean)
log_data$cv=100*log_data$sd/log_data$mean

order_data=log_data[order(log_data$cv,decreasing = T),]
ave_order_data=subset(order_data,order_data$mean>20)
top=ave_order_data[1:50,-c(22,23,24)]
plot=melt(top)
colnames(plot)=c("Motif","tissue","p.value")

lev=as.character(unique(plot$tissue))
order_lev=lev[order(lev,decreasing = T)]
plot$Tissue=factor(plot$tissue,levels = order_lev)

ggplot(plot,aes_(x = ~Motif, y = ~Tissue, size = ~p.value))+
    geom_point() +
    aes_string(color="p.value") +
    theme_bw()+
    theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
    scale_color_gradientn(colours = viridis(50,direction = -1,option = "D"))

##-------------------------------------------------------------##
##               Hclust of all TFs in samples                  ##
##-------------------------------------------------------------##
sub_data=subset(log_data,log_data$cv != "NaN")
tdata=t(sub_data[,-c(21:24)])
out.dist=dist(tdata,method="euclidean") 
out.hclust=hclust(out.dist,method="complete")
plot(out.hclust)   
