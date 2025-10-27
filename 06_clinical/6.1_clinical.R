#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

#install.packages("ggpubr")

library(limma)
library(ggpubr)
library(ComplexHeatmap)

expFile="rocSigExp.txt"       
cliFile="clinical.txt"       

# setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/06_clinical")

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

for(i in 1:length(rt[,1:ncol(rt)])){
    
    gene=colnames(rt)[i]
    
    tumorData=as.matrix(rt[gene])
    
    #rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
    
    data=avereps(tumorData)
    
    Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
    Type=factor(Type, levels=c("Low","High"))
    
    data1=cbind(as.data.frame(data), Type)
    data1=data1[order(data1[,gene]),] 
    #??ȡ?ٴ??????ļ?
    cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
    cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
    
    #?ϲ?????
    samSample=intersect(row.names(data), row.names(cli))
    data=data[samSample,,drop=F]
    cli=cli[samSample,,drop=F]
    rt1=cbind(data, cli)
    samSample2=intersect(row.names(data1), row.names(cli))
    data1=data1[samSample2,"Type",drop=F]
    cli=cli[samSample2,,drop=F]
    rt2=cbind(data1, cli)
    
    sigVec=c(gene)
    
    for(clinical in colnames(rt2[,2:ncol(rt2)])){
        data=rt2[c("Type", clinical)]
        colnames(data)=c("Type", "clinical")
        data=data[(data[,"clinical"]!="unknow"),]
        tableStat=table(data)
        stat=chisq.test(tableStat)
        pvalue=stat$p.value
        Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
        sigVec=c(sigVec, paste0(clinical, Sig))
    }
    
    colnames(rt2)=sigVec
    
    bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
             "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
             "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
             "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
             "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
             "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
    
    colorList=list()
    colorList[[gene]]=c("Low"="blue", "High"="red")
    j=0
    
    for(cli in colnames(rt2[,2:ncol(rt2)])){
        cliLength=length(levels(factor(rt2[,cli])))
        cliCol=bioCol[(j+1):(j+cliLength)]
        j=j+cliLength
        names(cliCol)=levels(factor(rt2[,cli]))
        cliCol["unknow"]="grey75"
        colorList[[cli]]=cliCol
    }
    
    #??????ͼ
    ha=HeatmapAnnotation(df=rt2, col=colorList)
    zero_row_mat=matrix(nrow=0, ncol=nrow(rt2))
    Hm=Heatmap(zero_row_mat, top_annotation=ha)
    
    #??????ͼ
    pdf(file=paste0(gene,"_heatmap.pdf"), width=7, height=5)
    draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
    #?ٴ??????Է?????????ͼ?ν???
    for(clinical in colnames(rt1[,2:ncol(rt1)])){
        data=rt1[c(gene, clinical)]
        colnames(data)=c(gene, "clinical")
        data=data[(data[,"clinical"]!="unknow"),]
        #???ñȽ???
        group=levels(factor(data$clinical))
        data$clinical=factor(data$clinical, levels=group)
        comp=combn(group,2)
        my_comparisons=list()
        for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
        #????????ͼ
        boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                          xlab=clinical,
                          ylab=paste(gene, " expression"),
                          legend.title=clinical)+ 
            stat_compare_means(comparisons = my_comparisons)
        #????ͼƬ
        pdf(file=paste0(gene,"_clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
        print(boxplot)
        dev.off()
    }
}


