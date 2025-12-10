---
title: "Differential analysis template (11/04/2016)"
author: "Zhen Gao"
date: "Mar 14, 2020"
output: 
    html_document: 
    fig_height: 10
fig_width: 10
---
    
    #### ZG updated based on 10082018 version
    #### Directory of output (modify code her                                                                                                                     e for different output directory or other projects)
    
    ```{r, echo=FALSE}
#rm(list=ls())
setwd("/media/code_12T/Dropbox/BBSR_P/2020/Karoline_Briegel")
Karoline_Briegel_dir = "/media/code_12T/Dropbox/BBSR_P/2020/Karoline_Briegel"
Date = Sys.Date()
Project = paste0("Karoline_Briegel_", Date)
count_dir = file.path(Karoline_Briegel_dir, "Count")
Out_dir = paste0(Karoline_Briegel_dir, "/Result/", Project, "/")
dir.create(paste0(Karoline_Briegel_dir, "/Result/"), showWarnings = F)
dir.create(Out_dir, showWarnings = F)
GSEA_folder = file.path(Karoline_Briegel_dir, "Result/GSEA/")
dir.create(GSEA_folder, showWarnings = F)
Out_dir
#ls()
```

#### Directory of sample information

```{r}
# data file
#cat(paste(dir(count_dir)[!grepl(".summary",dir(count_dir))], collapse="\n"))
SampleInfo_dir = Karoline_Briegel_dir
```

##### Packages 

```{r, message=FALSE, eval=TRUE}
# library("BiocInstaller")  # or source("https://bioconductor.org/biocLite.R") #'BiocInstaller' and 'biocLite()' are deprecated, use the 'BiocManager' CRAN package instead.
# library("BiocManager")         #### important changes
# BiocManager::install("edgeR")  #### important changes
# BiocManager::install("DESeq2")  #### important changes
# install.packages("devtools")
# install.packages("pheatmap")
# install.packages("plot3D")
# install.packages("edgeR")
# install.packages("BiocManager")
# install.packages("DESeq2")
require("devtools")
require("ggplot2")
require("pheatmap")
require("limma")
require(ggplot2)
require(grid)
require(plot3D)
require(multtest)
require(RColorBrewer)
require(survival)
require(limma)
require(edgeR)
require(plyr)
require(scales)
options(scipen=500)
require(DESeq2)
```

### Load count files and normalized data
#### Obtaining the sample information and cleaning and normalized raw data

```{r, eval=TRUE}
#---- Get count table----#
# Option 1: merge count
# merge raw count table
# merge_featureCount = function(DATA_dir, file_pattern, subtract_pattern = NULL){
#     File_name = dir(DATA_dir)
#     File_name = grep(file_pattern, File_name, value=T)
#     if(!is.null(subtract_pattern)){
#         File_name = File_name[!grepl(subtract_pattern, File_name)]
#     }
#     for(i in 1:length(File_name)){
#         #i=1
#         Temp_table = read.table(paste0(DATA_dir,"/",File_name[i]), header=T)
#         Temp_table = Temp_table[,c(1,7)]                                          #ZG, no $count in .count file
#         #dim(Temp_table)
#         #head(Temp_table)
#         if(i == 1){
#             Table_merged = Temp_table
#         }else{
#             Table_merged = merge(Table_merged, Temp_table, by = "Geneid", sort = F)
#         }
#     }
#     #sample_name = gsub(".*-(.*)_S.*","\\1", File_name)
#     #sample_name = gsub(".count", "", File_name)
#     rownames(Table_merged) = Table_merged[,"Geneid"]
#     Table_merged[,"Geneid"] = NULL
#     names(Table_merged) = File_name
#     Table_merged
# }
# rawdata.Karoline_Briegel = merge_featureCount(DATA_dir=count_dir,  file_pattern="STAR_", subtract_pattern = ".summary")
# rawdata.Karoline_Briegel_1   = read.table("/media/code_12T/Dropbox/BBSR_P/2020/Karoline_Briegel/readcount_genename.xls",      sep = "\t", header = T)  
# head(rawdata.Karoline_Briegel_1)
# dim(rawdata.Karoline_Briegel_1)
# head(rawdata.Karoline_Briegel_1)[, c(25:ncol(rawdata.Karoline_Briegel_1))]
# tail(rawdata.Karoline_Briegel_1)
# dim(rawdata.Karoline_Briegel_1)
########################################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# rawdata.Karoline_Briegel0   = read.csv("/media/code_12T/Dropbox/BBSR_P/2020/Karoline_Briegel/readcount_genename.csv",       header = T)  
# head(rawdata.Karoline_Briegel0)

rawdata.Karoline_Briegel   = read.csv("/media/code_12T/Dropbox/BBSR_P/2020/Karoline_Briegel/readcount_genename.csv",       header = T)  
head(rawdata.Karoline_Briegel)
dim(rawdata.Karoline_Briegel)
head(rawdata.Karoline_Briegel)[, c(25:ncol(rawdata.Karoline_Briegel))]
tail(rawdata.Karoline_Briegel)
dim(rawdata.Karoline_Briegel)

rownames(rawdata.Karoline_Briegel) = rawdata.Karoline_Briegel$geneID
rawdata.Karoline_Briegel = rawdata.Karoline_Briegel[,-c(1,26:ncol(rawdata.Karoline_Briegel))]
head(rawdata.Karoline_Briegel)

SampleInfo_dir=Karoline_Briegel_dir 
#Karoline_Briegel.Stype = read.csv(paste0(SampleInfo_dir,"/SampleID_INFO.csv"),header=T)  ##################### # the un-seen white space cause trouble
Karoline_Briegel.Stype = read.csv(paste0(SampleInfo_dir,"/SampleID_INFO-no-space.csv"),header=T)  #####################
Karoline_Briegel.Stype
#https://stackoverflow.com/questions/38874509/r-match-and-replace-column-names-by-data-frame
names(rawdata.Karoline_Briegel) <- Karoline_Briegel.Stype$New_ID[match(names(rawdata.Karoline_Briegel), Karoline_Briegel.Stype$ID)]

# Option 2: read merged count table prepared by featureCount_batch.0.1.R
# a. get final QC report
# please follow the code in ReadQC.batch.featureCounts.R 

# #---- Filter count----#
# # filtering genes with CPM 
rawdata.Karoline_Briegel.cpms <- cpm(rawdata.Karoline_Briegel)

class(rawdata.Karoline_Briegel.cpms)
head(rawdata.Karoline_Briegel.cpms)
dim(rawdata.Karoline_Briegel.cpms)
## filter out genes not expressed in at least two samples at 0.5 CPM
keep <- rowSums(rawdata.Karoline_Briegel.cpms  > 0.5)  >  2  ####################### 30 group ?? changed from 2, change it back to 2 06-122019
rawdata.Karoline_Briegel.filtrcpms <- rawdata.Karoline_Briegel[keep,]
head(rawdata.Karoline_Briegel.filtrcpms)
dim(rawdata.Karoline_Briegel.filtrcpms)

# sample ID information
Karoline_Briegel.Stype
#m.1 = match(colnames(rawdata.Karoline_Briegel), Karoline_Briegel.Stype$ID)
m.1 = match(colnames(rawdata.Karoline_Briegel), Karoline_Briegel.Stype$New_ID)
Karoline_Briegel.Stype = Karoline_Briegel.Stype[m.1, ]
Karoline_Briegel.Stype

# convert groups into factor
Karoline_Briegel.Stype$TYPE = factor(Karoline_Briegel.Stype$TYPE)

# class(rawdata.Karoline_Briegel.filtrcpms)
# dim(rawdata.Karoline_Briegel.filtrcpms)
# head(rawdata.Karoline_Briegel.filtrcpms)
# 
# rawdata.Karoline_Briegel.filtrcpms_colnames <- as.factor(colnames(rawdata.Karoline_Briegel.filtrcpms))
# rawdata.Karoline_Briegel.filtrcpms_colnames <- as.data.frame(rawdata.Karoline_Briegel.filtrcpms_colnames)
# 
# class(rawdata.Karoline_Briegel.filtrcpms_colnames)
# head(rawdata.Karoline_Briegel.filtrcpms_colnames)
# colnames(rawdata.Karoline_Briegel.filtrcpms_colnames) <- "New_ID"
# 
# 
# matched_Mtx = left_join(rawdata.Karoline_Briegel.filtrcpms_colnames, Karoline_Briegel.Stype, by = "New_ID")    # added 12-18-2018 merge may have order changes  #from dplyr
# colnames(rawdata.Karoline_Briegel.filtrcpms) <- matched_Mtx$SampleID
# head(rawdata.Karoline_Briegel.filtrcpms) 

write.csv(rawdata.Karoline_Briegel.filtrcpms, file.path(Out_dir, "rawdata.Karoline_Briegel.filtrcpms.csv"), row.names=T)    ###################

############# raw
# class(rawdata.Karoline_Briegel)
# dim(rawdata.Karoline_Briegel)
# head(rawdata.Karoline_Briegel)
# tail(rawdata.Karoline_Briegel)
# 
# rawdata.Karoline_Briegel_colnames <- as.factor(colnames(rawdata.Karoline_Briegel))
# rawdata.Karoline_Briegel_colnames <- as.data.frame(rawdata.Karoline_Briegel_colnames)
# 
# class(rawdata.Karoline_Briegel_colnames)
# head(rawdata.Karoline_Briegel_colnames)
# colnames(rawdata.Karoline_Briegel_colnames) <- "ID"
# matched_Mtx = left_join(rawdata.Karoline_Briegel_colnames, Karoline_Briegel.Stype, by = "ID")    # added 12-18-2018 merge may have order changes  #from dplyr
# colnames(rawdata.Karoline_Briegel) <- matched_Mtx$SampleID

write.csv(rawdata.Karoline_Briegel, file.path(Out_dir, "rawdata.Karoline_Briegel.csv"), row.names=T)                        #################

# data normalization
# note: Variance stabilizing apply to precleaning or after-cleaning (remove genes of all 0s) 
#       are the same
rawdata.Karoline_Briegel.wo0 = rawdata.Karoline_Briegel[rowSums(rawdata.Karoline_Briegel)>0, ]
head(rawdata.Karoline_Briegel)
head(rawdata.Karoline_Briegel.wo0)
 dim(rawdata.Karoline_Briegel.wo0)

 vsd.Karoline_Briegel.wo0     = varianceStabilizingTransformation(as.matrix(rawdata.Karoline_Briegel.wo0))
vsd.Karoline_Briegel         = varianceStabilizingTransformation(as.matrix(rawdata.Karoline_Briegel))
head(vsd.Karoline_Briegel.wo0)
head(vsd.Karoline_Briegel)
 dim(vsd.Karoline_Briegel)

# here we separate data into different batches for the differential analysis
type.notch = c("HCC1395_with_LBH_KD", "HCC1395_with_siCtrl", "MDA-MBA-MB-231_with_LBH_KD", "MDA-MBA-MB-231_with_231_siCtrl","MCF7_with_Vector_Control", "MCF7_with_LBH", "BT549_with_Vector_Control", "BT549_with_LBH")

rawdata.notch     = rawdata.Karoline_Briegel[, which(c(Karoline_Briegel.Stype$TYPE %in% type.notch))]
rawdata.notch.wo0 = rawdata.notch[rowSums(rawdata.notch)>0, ]
vsd.notch.wo0     = varianceStabilizingTransformation(as.matrix(rawdata.notch.wo0))
vsd.notch         = varianceStabilizingTransformation(as.matrix(rawdata.notch))
head(vsd.notch)

table_hs38ID = read.csv("/media/data_12T/Genome_Ref/gene_table/gene.table.GCv25.csv")  # local Ubuntu
#table_hs38ID = read.csv("/data/Genome_Ref/gene_table/gene_table_ucsc_hg38_04182017.csv")
head(table_hs38ID)
 dim(table_hs38ID)

table_hs38ID$ensembl_gene_id = sapply(strsplit(as.character(table_hs38ID$ensembl_gene_id), "\\."), "[[", 1 ) #or "[.]" https://stackoverflow.com/questions/26665100/how-to-use-the-strsplit-function-with-a-period

# mouse
#table_mm10ID = read.csv("/media/MyDATA/Genome_Ref/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/gene.table.GCvM12.csv", stringsAsFactors = F)
#table_mm10ID = read.csv("/data/Genome_Ref/gene_table/gene.table.GCvM12.csv", stringsAsFactors = F)
#head(table_mm10ID)

# output a normalized gene expression data
vsd.transformed.annot.data = merge(table_hs38ID, vsd.Karoline_Briegel, by.x = "ensembl_gene_id", by.y = "row.names", sort = F) # gene_id

 dim(vsd.transformed.annot.data)
head(vsd.transformed.annot.data)
write.csv(vsd.transformed.annot.data, file.path(Out_dir, "vsd.transformed.annot.data.csv"), row.names=T)

head(rawdata.Karoline_Briegel.filtrcpms)
 dim(rawdata.Karoline_Briegel.filtrcpms)

tmp_dds                 = DESeqDataSetFromMatrix(rawdata.Karoline_Briegel.filtrcpms, DataFrame(Karoline_Briegel.Stype), ~ 1)
tmp_dds                 = estimateSizeFactors(tmp_dds)
normSF.Karoline_Briegel = counts(tmp_dds, normalized=TRUE)
head(normSF.Karoline_Briegel)
 dim(normSF.Karoline_Briegel)

normSF.Karoline_Briegel.wAnno = merge(table_hs38ID, normSF.Karoline_Briegel, by.x = "ensembl_gene_id", by.y = "row.names", sort = F)
head(normSF.Karoline_Briegel.wAnno)
 dim(normSF.Karoline_Briegel.wAnno)
write.table(normSF.Karoline_Briegel.wAnno, file.path(Out_dir, "normSF.Karoline_Briegel.wAnno.xls"), row.names=F, sep = "\t")

#####################=========================== merge all raw and normalized reads as well as annotation ==================================================================================================================
head(normSF.Karoline_Briegel.wAnno)
 dim(normSF.Karoline_Briegel.wAnno)
annoted_raw_and_normalized_reads <- merge(normSF.Karoline_Briegel.wAnno, rawdata.Karoline_Briegel,by.x= "ensembl_gene_id",  by.y= 0)   # "row.names" or the number 0 specifies the row names
head(annoted_raw_and_normalized_reads)
 dim(annoted_raw_and_normalized_reads)

write.table(annoted_raw_and_normalized_reads, file.path(Out_dir,"annoted_raw_and_normalized_reads.xls"), sep = "\t", quote=F, row.names = FALSE) # ".xls" "csv" "\t", the name to be true <==03-13row.names= T, 
getwd()
# interesting_genes <- c("IL2ra", "FoxP3")
# head(rawdata.Karoline_Briegel.filtrcpms)
#  dim(rawdata.Karoline_Briegel.filtrcpms)
# annoted_raw_and_normalized_reads <- merge(normSF.Karoline_Briegel.wAnno, rawdata.Karoline_Briegel,by.x= "ensembl_gene_id",  by.y= 0)   # "row.names" or the number 0 specifies the row names 

```

# clustring with the 500 most variable genes

```{r}
# myvars <- apply(normSF.Karoline_Briegel, 1, var, na.rm=TRUE) 
# myvars <- sort(myvars, decreasing=TRUE) 
# myvars <- myvars[1:500] 
# head(myvars)
# length(myvars)
# normSF.Karoline_Briegel_500 <- normSF.Karoline_Briegel[names(myvars), ] 
#  dim(normSF.Karoline_Briegel_500) 
# head(normSF.Karoline_Briegel_500) 
# # Create a matrix
# hclust_matrix <- normSF.Karoline_Briegel_500 %>%  as.matrix()
# head(hclust_matrix)
# hclust_matrix <- hclust_matrix %>% 
#   # transpose the matrix so genes are as columns
#   t() %>% 
#   # apply scalling to each column of the matrix (genes)
#   scale() 
# 
# sample_dist   <- dist(hclust_matrix)
# sample_hclust <- hclust(sample_dist, method = "complete")
# 
# plot(sample_hclust)
# abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
```


clustring based on all 20K genes

```{r}
# https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/
# install.packages("factoextra")
require("factoextra")
hclust_matrix <- normSF.Karoline_Briegel %>%  as.matrix()
head(hclust_matrix)
 dim(hclust_matrix)

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() 

sample_dist   <- dist(hclust_matrix)
sample_hclust <- hclust(sample_dist, method = "complete")

# plot(sample_hclust)
# abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
fviz_dend(sample_hclust, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
          )


dev.copy(pdf, file=paste0(Out_dir,"/","sample_clustering.pdf"), width=8, height=8)
dev.off() 



```

#### Make heatmap and PCA plot

```{r, eval=FALSE}
##################

# heatmap and pca
# PCA and heatmap
PCA_3D = function(Data, output_pca, out_dir, g1_level = NULL, g2_level = NULL, legend.add=F){
    #clean data column name
    colnames(Data) = gsub("STAR_", "", colnames(Data))
    colnames(Data) = gsub(".count", "", colnames(Data))
    colnames(Data) = gsub("_R1", "", colnames(Data))
    
    hmcol<-rev(colorRampPalette(brewer.pal(9, "Set1"))(256))
    pch = 16
    if(is.null(g1_level)){
        type_level = pch
        col_level = "black"
    }else{
        TEMP = factor(g1_level)
        uniq_label_1 =  levels(TEMP)
        levels(TEMP) = hmcol[ceiling(seq(length.out=length(levels(TEMP)),from=1,to=256))]
        col_level = as.character(TEMP)
        uniq_col = levels(TEMP)
        if(!is.null(g2_level)){
            TEMP = factor(g2_level)
            uniq_label_2 =  levels(TEMP)
            levels(TEMP) = as.character(c(19,17,15,18,3,4,7,8,9,10)[1:length(levels(TEMP))])
            type_level = as.numeric(as.character(TEMP))
            uniq_type = as.numeric(levels(TEMP))
        }else{
            type_level = pch
        }
    }
    
    #png(file=paste0(out_dir,"/",output_pca), width=2000, height=2000, res = 300)
    fig.w = 7
    fig.h = 7
    par(mar=c(6.1, 4.1, 1.5, 3.1), xpd=TRUE)
    par(oma=c(0,0,0,0))
    Data.pca = prcomp(t(Data), center = T, scale.=F)
    with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, 
                                           type="h", lty.hplot=2, cex=1,
                                           ticktype = "detailed", nticks = 2,
                                           bty = "b2", #bty="b2", 
                                           xlab="PC 1",	ylab="PC 2",zlab="PC 3", 
                                           theta = -50, phi = 40, 
                                           pch=type_level,
                                           col=col_level,
                                           main = "Principal component analysis")
    )
    if(!is.null(g1_level)){
        legend("right", legend = uniq_label_1, pch=pch, 
               col = uniq_col, bty= "n", inset=c(0.2,0.3), xjust=0,
               cex=0.8)
    }
    if(!is.null(g2_level)){
        legend("bottomright", legend = uniq_label_2, pch=uniq_type, 
               col = "grey", bty= "n", inset=c(0.2,-0.1), xjust=0,
               cex=0.8)
    }
    if(!legend.add){
        sample_label = colnames(Data)
    }else{
        sample_label = 1:ncol(Data)
    }
    with(data.frame(Data.pca$x), text3D(x=PC1, y=PC2, 
                                        z=PC3, labels = integer(length(PC1)), col = "black", add=TRUE, colkey = FALSE, adj=1.7, cex=0.5) #sample_label ============== remove label
    )
    if(legend.add){
        legend("bottom", legend = paste0(1:ncol(Data), ": ", colnames(Data)), pch=0, yjust = 0, 
               col = "white", cex=0.3, bty= "n", ncol = 2, inset=c(-0.1,-0.2))
        fig.h = 10
    }
    dev.copy(pdf, file=paste0(out_dir,"/",output_pca), width=fig.w, height=fig.h)
    dev.off()
}
PCA_3D(vsd.Karoline_Briegel.wo0, "PCA_Karoline_Briegel3D.pdf", Out_dir, Karoline_Briegel.Stype$TYPE)   # note ,should be pdf

# pca for selected variables with 2 ID 

PCA_2D = function(Data, output_pca, out_dir, g_level, g2_level=NULL, clustering = NULL, label = FALSE,  label_info = NA, X_lim = NULL, Y_lim = NULL, g_lab = "", g2_lab = NA){
    
    tmp.pca = prcomp(t(Data), center = T, scale.=F)
    # tmp_pca_data_frame = as.data.frame(tmp.pca$x)  
    # tmp_pca_data_frame$ID = row.names(tmp_pca_data_frame)  
    # selected_YBZhang.Stype = YBZhang.Stype[, c("ID", "SampleID")] 
    # output_PCA = left_join(selected_YBZhang.Stype, tmp_pca_data_frame, by = "ID")       ################ added 12-12-2018 merge may have order changes  #from dplyr
    # write.csv(output_PCA[, c(1:4)], paste0(Out_dir, output_pca, ".csv"))
    
    Tmp_1 = summary(tmp.pca)
    var.pc1 = percent(Tmp_1$importance[2, "PC1"])
    var.pc2 = percent(Tmp_1$importance[2, "PC2"])
    
    TEMP = factor(as.character(g_level))
    uniq_label =  levels(TEMP)
    levels(TEMP) = rev(colorRampPalette(brewer.pal(9, "Set1"))(256))[ceiling(seq(length.out=length(uniq_label),from=1,to=256))]
    
    TEMP.PCA = data.frame(tmp.pca$x, T_level = factor(g_level))

    if(label){ 
    Temp.plot = ggplot(data = TEMP.PCA, aes(x = PC1, y = PC2, label = label_info))+  # HERE THE LABEL IS NECESSARY FOR ggrepel by ZG
            geom_text_repel( size=1.5) # point.padding = NA, segment.size = 0.2, size=1.5
    } else {
        Temp.plot = ggplot(data = TEMP.PCA, aes(x = PC1, y = PC2))
    }

    GG = Temp.plot + 
        coord_fixed() +  
        scale_color_manual(name= g_lab,values = levels(TEMP)) + 
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(colour = "grey", fill=NA, size=1), 
              legend.position = "right", legend.direction = "vertical", legend.key = element_blank()) +
        xlab(paste0("Principal component 1 (", var.pc1, ")")) + 
        ylab(paste0("Principal component 2 (", var.pc2, ")")) 
    
    if(!is.null(g2_level)){
        TEMP2 = factor(g2_level)
        levels(TEMP2) = as.numeric(c(19,17,15,18,3,4,7,8,9,10)[1:length(levels(TEMP2))])
        GG = GG + geom_point(aes(colour = T_level, shape = g2_level)) + scale_shape_manual(name=g2_lab, values = as.numeric(levels(TEMP2)))
    }else{
        GG = GG + geom_point(aes(colour = T_level))
    }
    # if(label){
    #     GG = GG + geom_text(aes(label=label), size=1.5, hjust = -0.2, nudge_x = 1)   
    # }
    if(!is.null(clustering)){
        GG = GG + stat_ellipse(aes(group = Cluster), color="black", linetype = 2, type = "t")
    }
    if(!is.null(X_lim)){
        GG = GG + coord_cartesian(xlim = X_lim)
    }
    
    output_pca = paste0(output_pca, ".pdf")   ### added 
    ggsave(file.path(out_dir, output_pca),
           GG,
           device = "pdf",
           width = 8,   # 6
           height = 6)  #4.25
}

PCA_2D(Data=vsd.Karoline_Briegel.wo0, 
       output_pca="PCA2D_Karoline_Briegel_all_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# Karoline_Briegel.Stype$New_ID
# overexpressed_ID = c("B_C1", "B_C2", " B_C3","B_LBH1","B_LBH2" ,"B_LBH3", "H_C1","H_C2"," H_C3","H_KD1", "H_KD2","H_KD3","MD231_C1","MD231_C2", "MD231_C3","MD231_KD1","MD231_KD2","MD231_KD3", "ME_C1","ME_C2","ME_C3", "ME_LBH1","ME_LBH2"," ME_LBH3")

HCC_and_KD_ID     = c("H_C1",     "H_C2",     "H_C3",     "H_KD1",     "H_KD2",     "H_KD3"    )
MDA_and_KD_ID     = c("MD231_C1", "MD231_C2", "MD231_C3", "MD231_KD1", "MD231_KD2", "MD231_KD3")
combined_KD_ID    = c("H_C1",     "H_C2",     "H_C3",     "H_KD1",     "H_KD2",     "H_KD3",   "MD231_C1", "MD231_C2", "MD231_C3", "MD231_KD1", "MD231_KD2", "MD231_KD3")

BT_and_LBH_OV_ID  = c("B_C1",     "B_C2",     "B_C3",     "B_LBH1",    "B_LBH2" ,   "B_LBH3"   )
MCF_and_LBH_OV_ID = c("ME_C1",    "ME_C2",    "ME_C3",    "ME_LBH1",   "ME_LBH2",   "ME_LBH3"  )                             
combined_OV_ID    = c("ME_C1",    "ME_C2",    "ME_C3",    "ME_LBH1",   "ME_LBH2",   "ME_LBH3", "B_C1",     "B_C2",     "B_C3",     "B_LBH1",    "B_LBH2" ,   "B_LBH3")

#head(vsd.Karoline_Briegel.wo0)       
# 1) HCC_and_KD
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% HCC_and_KD_ID)], 
       output_pca="PCA2D_Karoline_Briegel_HCC_and_KD_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% HCC_and_KD_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% HCC_and_KD_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# 2) MDA_and_KD
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% MDA_and_KD_ID)], 
       output_pca="PCA2D_Karoline_Briegel_MDA_and_KD_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% MDA_and_KD_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% MDA_and_KD_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# 3) combined_KD
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% combined_KD_ID)], 
       output_pca="PCA2D_Karoline_Briegel_combined_KD_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% combined_KD_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% combined_KD_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# 4) MCF_and_LBH_OV
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% MCF_and_LBH_OV_ID)], 
       output_pca="PCA2D_Karoline_Briegel_MCF_and_LBH_OV_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% MCF_and_LBH_OV_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% MCF_and_LBH_OV_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# 5) BT_and_LBH_OV
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% BT_and_LBH_OV_ID)], 
       output_pca="PCA2D_Karoline_Briegel_BT_and_LBH_OV_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% BT_and_LBH_OV_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% BT_and_LBH_OV_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf

# combined_OV
PCA_2D(Data=vsd.Karoline_Briegel.wo0[,which(colnames(vsd.Karoline_Briegel.wo0) %in% combined_OV_ID)], 
       output_pca="PCA2D_Karoline_Briegel_combined_OV_groups", 
       out_dir= Out_dir, 
       g_level=Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% combined_OV_ID), ]$TYPE,   
       label = T,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       label_info = Karoline_Briegel.Stype[which(Karoline_Briegel.Stype$New_ID %in% combined_OV_ID), ]$New_ID,     #  point txt label  label = Karoline_Briegel.Stype$SampleID,  
       g_lab = "Treatment-Group")  # .pdf


# after check the PCA ME_C1 (MCF7 control1) is outliar 
# after check the PCA MD231_KD3 might be a outliar
# Decided to delete ME_C1 (MCF7 control1) from the analysis 
# after clustring the H_C3 was removed as outlier
```



### Run DESeq analysis

```{r, eval=TRUE}
# run DESeq analysis
# DESeq_2group: comparing trt with control
# DESeq_mlutigroupLRT: overall test whether there is difference among groups
DESeq_2group = function(DATA_raw, Sample_ty1, Sample_ty2, Sample_info){
    cat("Updated 04/08/2017 *use pooled variance from ALL samples*
        Note: Universal gene name must be gene symbol with column identifying name 'gene_id'\n")
    cat("   Sample_info must be supplied with the column format 'ID' corresponds colnames(DATA_raw)
        and 'TYPE' corresponds Sample_ty1 and Sample_ty2\n")
    options(scipen=500)
    cat("Matching Sample ID...\n")
    w.1 = match(Sample_info$ID, colnames(DATA_raw))
    cat("  all selected sample names:\n")
    print(colnames(DATA_raw)[w.1])
    Data_s = DATA_raw[, w.1]
    if(length(grep("TYPE1", colnames(Sample_info)))==0){
        Sample_id1 = which(Sample_info$TYPE==Sample_ty1)
        Sample_id2 = which(Sample_info$TYPE==Sample_ty2)
        condition = factor(Sample_info$TYPE)
        dds = DESeqDataSetFromMatrix(Data_s, DataFrame(condition), ~ condition)
    }else{
        cat("Samples stratified into multiple types using TYPE2 information\n")
        Sample_id1 = which(Sample_info$TYPE1==Sample_ty1)
        Sample_id2 = which(Sample_info$TYPE1==Sample_ty2)
        condition = factor(Sample_info$TYPE1)
        type2 = factor(Sample_info$TYPE2)
        dds = DESeqDataSetFromMatrix(Data_s, DataFrame(data.frame(condition, type2)), ~ condition + type2)
    }
    return(results(DESeq(dds), c("condition", Sample_ty1, Sample_ty2)))
}

# automatic combine differential results with annotated and normalized data
DESeq_2group_wAnno = function(DATA_raw, Sample_ty1, Sample_ty2, Sample_info, Anno_table, out_dir, 
                              coldrop_DE = c("lfcSE"), sort.by = "padj", use_log2 = T){
    cat("Note: Updated 04/08/2017 *use pooled variance from ALL samples*
        Universal gene name must be gene symbol with column identifying name 'gene_id'\n")
    cat("   Sample_info must be supplied with the column format 'ID' corresponds colnames(DATA_raw)
        and 'TYPE' corresponds Sample_ty1 and Sample_ty2\n")
    options(scipen=500)
    cat("Matching Sample ID...\n")
    w.1 = match(Sample_info$ID, colnames(DATA_raw))
    cat("  all selected sample names:\n")
    
    print(colnames(DATA_raw)[w.1])
    Data_s = DATA_raw[, w.1]
    if(length(grep("TYPE1", colnames(Sample_info)))==0){
        Sample_id1 = which(Sample_info$TYPE==Sample_ty1)
        Sample_id2 = which(Sample_info$TYPE==Sample_ty2)
        condition = factor(Sample_info$TYPE)
        dds = DESeqDataSetFromMatrix(Data_s, DataFrame(condition), ~ condition)
    }else{
        cat("Samples stratified into multiple types using TYPE2 information\n")
        Sample_id1 = which(Sample_info$TYPE1==Sample_ty1)
        Sample_id2 = which(Sample_info$TYPE1==Sample_ty2)
        condition = factor(Sample_info$TYPE1)
        type2 = factor(Sample_info$TYPE2)
        dds = DESeqDataSetFromMatrix(Data_s, DataFrame(data.frame(condition, type2)), ~ condition + type2)
    }
    Temp_DE = results(DESeq(dds), c("condition", Sample_ty1, Sample_ty2))
    if(!use_log2){
        colnames(Temp_DE)[grep("log2FoldChange", colnames(Temp_DE))] = "FoldChange"
        Temp_DE$FoldChange = 2^Temp_DE$FoldChange
        Temp_DE$FoldChange = ifelse(Temp_DE$FoldChange>=1, Temp_DE$FoldChange, -1/Temp_DE$FoldChange)
    }else{
        colID_log2FC = grep("log2FoldChange", colnames(Temp_DE))
        Temp_DE = data.frame(Temp_DE[, 1:colID_log2FC], FoldChange = 2^Temp_DE$log2FoldChange, Temp_DE[, (colID_log2FC+1):ncol(Temp_DE)], check.names = F)
    }
    
    Temp_DE = data.frame(gene_id_DE = rownames(Temp_DE), Temp_DE[, !(colnames(Temp_DE) %in% coldrop_DE)], check.names = F)
    
    Sample.name.match = match(c(as.character(Sample_info$ID[Sample_id1]), as.character(Sample_info$ID[Sample_id2])), colnames(DATA_raw))
    Temp_DEwCount = merge(Temp_DE, data.frame(gene_id_DE = rownames(DATA_raw), DATA_raw[, Sample.name.match], check.names = F), by="gene_id_DE", all.x=F, all.y=T, sort=F)
    
    # determine which gene annotation files to use
    if(colnames(Anno_table)[1]=="ensembl_gene_id"){
        Temp_DEwAnno = merge(Anno_table, Temp_DEwCount, by.x = "ensembl_gene_id", by.y ="gene_id_DE", all.x=F, all.y = T, sort=F)
    }else{
        Temp_DEwAnno = merge(Anno_table, Temp_DEwCount, by.x = "gene_id", by.y="gene_id_DE", all.x=F, all.y = T, sort=F)
    }
    
    #print(colnames(Temp_DEwAnno))
    # sort data
    if(sort.by=="padj"){
        write.table(Temp_DEwAnno[order(Temp_DEwAnno$padj, Temp_DEwAnno$pvalue),], paste0(out_dir, "/", Sample_ty1,"vs",Sample_ty2,".csv"), row.names=F, sep = ",", quote=F) # change to csv 03-14-2020
    }else{
        write.table(Temp_DEwAnno, paste0(out_dir, "/", Sample_ty1,"vs",Sample_ty2,".csv"),  row.names=F, sep = ",", quote=F) # ".xls"),  row.names=F, sep = "\t" 03-14-2020
    }
}

BT_and_LBH_OV_ID        = c("B_C1",     "B_C2",     "B_C3",     "B_LBH1",    "B_LBH2" ,   "B_LBH3"   )
HCC_and_KD_legal_ID     = c("H_C1",     "H_C2",                 "H_KD1",     "H_KD2",     "H_KD3"    ) #"H_C3", 

MCF_and_LBH_OV_legal_ID = c(            "ME_C2",    "ME_C3",    "ME_LBH1",   "ME_LBH2",   "ME_LBH3"  )  
MDA_and_KD_ID           = c("MD231_C1", "MD231_C2", "MD231_C3", "MD231_KD1", "MD231_KD2", "MD231_KD3")
                           


legal_ID    = c(HCC_and_KD_legal_ID, MDA_and_KD_ID, MCF_and_LBH_OV_legal_ID, BT_and_LBH_OV_ID ) 
                
rawdata.Karoline_Briegel.filtrcpms_with_legal_ID = rawdata.Karoline_Briegel.filtrcpms[,which(colnames(rawdata.Karoline_Briegel.filtrcpms) %in% legal_ID)]
 dim(rawdata.Karoline_Briegel.filtrcpms)
head(rawdata.Karoline_Briegel.filtrcpms)
head(rawdata.Karoline_Briegel.filtrcpms_with_legal_ID)
 dim(rawdata.Karoline_Briegel.filtrcpms_with_legal_ID)

### prepare the table 
#rawdata.Karoline_Briegel.filtrcpms_with_HCC_and_KD_legal_ID = rawdata.Karoline_Briegel.filtrcpms[,which(colnames(rawdata.Karoline_Briegel.filtrcpms) %in% HCC_and_KD_legal_ID)]
Karoline_Briegel.Stype1     = Karoline_Briegel.Stype                               
Karoline_Briegel.Stype1$ID  = NULL
names(Karoline_Briegel.Stype1)[names(Karoline_Briegel.Stype1) == "New_ID"] <- "ID"
Karoline_Briegel.Stype1          ##########################################======================== make the ID consistent accros samples ME_L1  ==>  ME_LBH1 to indicated oveexoressuion LDH or MD231L_1 ==> MD231L_KD to indicate knowkdown

head(rawdata.Karoline_Briegel.filtrcpms_with_legal_ID)
 dim(rawdata.Karoline_Briegel.filtrcpms_with_legal_ID)

# 1) HCC_and_KD_legal_ID
Karoline_Briegel.Stype2     = Karoline_Briegel.Stype1[which(Karoline_Briegel.Stype1$ID %in% HCC_and_KD_legal_ID), ]
DESeq_2group_wAnno(DATA_raw = rawdata.Karoline_Briegel.filtrcpms_with_legal_ID, Sample_ty1="HCC1395_with_LBH_KD", Sample_ty2="HCC1395_with_siCtrl",           Sample_info=Karoline_Briegel.Stype2, Anno_table=table_hs38ID, out_dir=Out_dir)
# 2) MDA_and_KD_ID
Karoline_Briegel.Stype3     = Karoline_Briegel.Stype1[which(Karoline_Briegel.Stype1$ID %in% MDA_and_KD_ID), ]
DESeq_2group_wAnno(DATA_raw = rawdata.Karoline_Briegel.filtrcpms_with_legal_ID, Sample_ty1="MDA-MBA-MB-231_with_LBH_KD", Sample_ty2="MDA-MBA-MB-231_with_231_siCtrl",Sample_info=Karoline_Briegel.Stype3, Anno_table=table_hs38ID, out_dir=Out_dir)
# 3) MCF_and_LBH_OV_legal_ID
Karoline_Briegel.Stype4     = Karoline_Briegel.Stype1[which(Karoline_Briegel.Stype1$ID %in% MCF_and_LBH_OV_legal_ID), ]
DESeq_2group_wAnno(DATA_raw = rawdata.Karoline_Briegel.filtrcpms_with_legal_ID, Sample_ty1="MCF7_with_LBH",       Sample_ty2="MCF7_with_Vector_Control",   Sample_info=Karoline_Briegel.Stype4, Anno_table=table_hs38ID, out_dir=Out_dir)
# 4) BT_and_LBH_OV_ID
Karoline_Briegel.Stype5     = Karoline_Briegel.Stype1[which(Karoline_Briegel.Stype1$ID %in% BT_and_LBH_OV_ID), ]
DESeq_2group_wAnno(DATA_raw = rawdata.Karoline_Briegel.filtrcpms_with_legal_ID, Sample_ty1="BT549_with_LBH",      Sample_ty2="BT549_with_Vector_Control",    Sample_info=Karoline_Briegel.Stype5, Anno_table=table_hs38ID, out_dir=Out_dir)

```

# heat maps
```{r, eval=FALSE}
# two group heatmap comparison, 
#  select top DE 100 genes
#  and separate them into up and down groups
heatmap_2group = function(Data, Pvalue, log2FC, output_heatmap, out_dir, cutoffs = NULL, 
                          padj_thrshd = NULL, high_id = NULL, low_id = NULL, dist_fun = "euclidean",
                          gene_order_by = c("cluster", "fc", "none"),
                          col_range=c("redblue", "redgreen"), main=NULL, return_matrix=F){
    require(gplots)
    col_range = col_range[1]
    if(col_range == "redgreen"){
        hmcol<-rev(redgreen(50))[-seq(35, 35)] 
    }else if(col_range == "redblue"){
        hmcol<-rev(redblue(50))[-seq(35, 35)]
    }
    message("If use cutoffs, the format is cutoffs=c(pval_cutoff, FC_cutoff)\n
            If use high and low, specify the ids of columns\n")
    #hmcol<-rev(colorRampPalette(brewer.pal(10, "RdYlGn"))(256))
    #hmcol<-rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
    if(is.null(high_id)&is.null(low_id)){
        Data = Data
        heatmap_colorder = 1:ncol(Data)
    }else{
        Data = Data[, c(high_id, low_id)]
        heatmap_colorder = 1:ncol(Data)
    }
    if(!is.null(padj_thrshd)){
        DE.sort = which(Pvalue < padj_thrshd)
        title_lab = paste0("Features with significant DE for ", output_heatmap, " (adj.p < ",  padj_thrshd ,")")
    }else if(is.null(cutoffs)){
        DE.sort = order(Pvalue)[1:100]
        title_lab = paste0("Top 100 most DE features\n(", output_heatmap, ")")
    }else{
        DE.sort = which((Pvalue < cutoffs[1])&(abs(log2FC) > cutoffs[2]))
        title_lab = paste0("Features with significant DE\n(", output_heatmap, ")")
    }
    if(!is.null(main)){
        title_lab = main
    }
    if(length(DE.sort)<=100){
        labelrow = NULL
    }else{
        labelrow = NA
    }
    TEMP.data = data.frame(Data[DE.sort,])
    TEMP.up = TEMP.data[log2FC[DE.sort]>=0,]
    TEMP.dn = TEMP.data[log2FC[DE.sort]<0,]
    gene_order_by = gene_order_by[1]
    
    if(gene_order_by=="cluster"){ # default
        if(nrow(TEMP.up)>2){
            Data.up.hc = hclust(dist(TEMP.up), method="complete")  
            rowInd.up <- order.dendrogram(as.dendrogram(Data.up.hc))
            Temp1 = TEMP.data[log2FC[DE.sort]>=0,][rowInd.up,]
        }else{
            Temp1 = TEMP.data[log2FC[DE.sort]>=0,]
        }
        if(nrow(TEMP.dn)>2){
            Data.dn.hc = hclust(dist(TEMP.dn), method="complete")
            rowInd.dn <- order.dendrogram(as.dendrogram(Data.dn.hc))    
            Temp2 = TEMP.data[log2FC[DE.sort]<0,][rowInd.dn,]
        }else{
            Temp2 = TEMP.data[log2FC[DE.sort]<0,]
        }
        Temp = rbind(Temp1, Temp2)
    }else if(gene_order_by=="fc"){
        o.1 = order(log2FC[DE.sort], decreasing = T)
        Temp = TEMP.data[o.1,]
    }else if(gene_order_by=="none"){
        Temp = data.frame(Data)
    }
    
    Data.sorted = as.matrix(Temp)
    
    #Data.sorted = t(apply(Data.sorted, 1, scale)) # add scaling by row for clustering????????
    
    rownames(Data.sorted) = rownames(Temp)
    colnames(Data.sorted) = colnames(Data)
    
    #write.table(data.frame(rownames(Data.sorted)), file = paste0(out_dir,"/",output_heatmap,".genelist.csv"), row.names=F, col.names = F, quote=F)
    if(nrow(Data.sorted)>1){
        pdf(file=paste0(out_dir,"/",output_heatmap,".pdf"), width=2000/400, height=ceiling(3600/400*(0.66+nrow(Data.sorted)/300)))
        par(cex.main=0.5)
        heatmap.2(Data.sorted, hclustfun = function(h)hclust(h, method="complete"), col=hmcol, 
                  distfun = function(d)dist(d, method=dist_fun), 
                  Colv = heatmap_colorder, 
                  reorderfun = function(d, w)reorder(d, w, mean),
                  Rowv = F, dendrogram = "column", srtCol = 45, labRow = labelrow, cexRow = 0.5,
                  scale="row",key=TRUE, keysize=0.5, symkey=FALSE,density.info="none", trace="none",cexCol=0.6,
                  margins=c(4,5), lwid = c(1,4), lhei = c(0.9,5.8), key.par=list(mar=c(4, 1, 2, 0)),
                  main = title_lab) #
        dev.off()
    }
    
    if(return_matrix){
        if(return_matrix){
            Data.sorted
        }
    }
}

# DE_dir = paste0(Out_dir, "Differential_Genes")
BT549_with_LBHvsBT549_with_Vector_Control  = read.csv(file.path(Out_dir,"BT549_with_LBHvsBT549_with_Vector_Control.csv"),                 header = T) ####$modified by ZG###  read.table, xls sep = "\t"
HCC1395_with_LBH_KDvsHCC1395_with_siCtrl   = read.csv(file.path(Out_dir,"HCC1395_with_LBH_KDvsHCC1395_with_siCtrl.csv"),                  header = T) ####$modified by ZG### sep = "\t", 
MCF7_with_LBHvsMCF7_with_Vector_Control    = read.csv(file.path(Out_dir,"MCF7_with_LBHvsMCF7_with_Vector_Control.csv"),                   header = T) ####$modified by ZG###
MDA231_with_LBH_KDvsMDA231_with_231_siCtrl = read.csv(file.path(Out_dir,"MDA-MBA-MB-231_with_LBH_KDvsMDA-MBA-MB-231_with_231_siCtrl.csv"),header = T) ####$modified by ZG###

head(BT549_with_LBHvsBT549_with_Vector_Control)
 dim(BT549_with_LBHvsBT549_with_Vector_Control)

colnames(vsd.transformed.annot.data)
    head(vsd.transformed.annot.data)
     dim(vsd.transformed.annot.data)

legal_colnames =c("ensembl_gene_id" ,"gene_id", legal_ID)
vsd_BT549_with_LBHvsBT549_with_Vector_Control  = left_join(BT549_with_LBHvsBT549_with_Vector_Control[ ,  c("ensembl_gene_id","log2FoldChange","pvalue")] , vsd.transformed.annot.data[, legal_colnames] )
vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl   = left_join(HCC1395_with_LBH_KDvsHCC1395_with_siCtrl[ ,   c("ensembl_gene_id","log2FoldChange","pvalue")] , vsd.transformed.annot.data[, legal_colnames] )
vsd_MCF7_with_LBHvsMCF7_with_Vector_Control    = left_join(MCF7_with_LBHvsMCF7_with_Vector_Control[ ,    c("ensembl_gene_id","log2FoldChange","pvalue")] , vsd.transformed.annot.data[, legal_colnames] )
vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl = left_join(MDA231_with_LBH_KDvsMDA231_with_231_siCtrl[ , c("ensembl_gene_id","log2FoldChange","pvalue")] , vsd.transformed.annot.data[, legal_colnames] )

head(vsd_BT549_with_LBHvsBT549_with_Vector_Control)
 dim(vsd_BT549_with_LBHvsBT549_with_Vector_Control)

vsd_BT549_with_LBHvsBT549_with_Vector_Control  = vsd_BT549_with_LBHvsBT549_with_Vector_Control[ which( vsd_BT549_with_LBHvsBT549_with_Vector_Control$gene_id!="NA"), ]  
vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl   = vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl[  which(  vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl$gene_id!="NA"), ]
vsd_MCF7_with_LBHvsMCF7_with_Vector_Control    = vsd_MCF7_with_LBHvsMCF7_with_Vector_Control[   which(   vsd_MCF7_with_LBHvsMCF7_with_Vector_Control$gene_id!="NA"), ]
vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl = vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl[which(vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl$gene_id!="NA"), ]

head(vsd_MCF7_with_LBHvsMCF7_with_Vector_Control )
 dim(vsd_MCF7_with_LBHvsMCF7_with_Vector_Control )

rownames(vsd_BT549_with_LBHvsBT549_with_Vector_Control)  <- make.names( vsd_BT549_with_LBHvsBT549_with_Vector_Control$gene_id, unique=TRUE)
rownames(vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl)   <- make.names(  vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl$gene_id, unique=TRUE)
rownames(vsd_MCF7_with_LBHvsMCF7_with_Vector_Control)    <- make.names(   vsd_MCF7_with_LBHvsMCF7_with_Vector_Control$gene_id, unique=TRUE)
rownames(vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl) <- make.names(vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl$gene_id, unique=TRUE)

head(vsd_BT549_with_LBHvsBT549_with_Vector_Control)  
head(vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl)   
head(vsd_MCF7_with_LBHvsMCF7_with_Vector_Control) 
head(vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl) 


BT_and_LBH_OV_ID        = c("B_C1",     "B_C2",     "B_C3",     "B_LBH1",    "B_LBH2" ,   "B_LBH3"   )
HCC_and_KD_legal_ID     = c("H_C1",     "H_C2",                 "H_KD1",     "H_KD2",     "H_KD3"    ) #"H_C3", 
MCF_and_LBH_OV_legal_ID = c(            "ME_C2",    "ME_C3",    "ME_LBH1",   "ME_LBH2",   "ME_LBH3"  )  
MDA_and_KD_ID           = c("MD231_C1", "MD231_C2", "MD231_C3", "MD231_KD1", "MD231_KD2", "MD231_KD3")

heatmap_2group(Data=vsd_BT549_with_LBHvsBT549_with_Vector_Control[ ,BT_and_LBH_OV_ID],
               Pvalue=vsd_BT549_with_LBHvsBT549_with_Vector_Control$pvalue,
               log2FC=vsd_BT549_with_LBHvsBT549_with_Vector_Control$log2FoldChange,
               output_heatmap="BT549_with_LBH_vs_BT549_with_Vector_Control",
               out_dir=Out_dir)

heatmap_2group( Data=vsd_MCF7_with_LBHvsMCF7_with_Vector_Control[ ,MCF_and_LBH_OV_legal_ID], 
                Pvalue=MCF7_with_LBHvsMCF7_with_Vector_Control$pvalue, 
                log2FC=MCF7_with_LBHvsMCF7_with_Vector_Control$log2FoldChange, 
                output_heatmap="MCF7_with_LBHvsMCF7_with_Vector_Control", 
                out_dir=Out_dir)

heatmap_2group(Data   = vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl[ ,HCC_and_KD_legal_ID],
               Pvalue = vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl$pvalue,
               log2FC = vsd_HCC1395_with_LBH_KDvsHCC1395_with_siCtrl$log2FoldChange,
               output_heatmap="HCC1395_with_LBH_KDvsHCC1395_with_siCtrl",
               out_dir=Out_dir)


heatmap_2group(Data   = vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl[ ,MDA_and_KD_ID],
               Pvalue = vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl$pvalue,
               log2FC = vsd_MDA231_with_LBH_KDvsMDA231_with_231_siCtrl$log2FoldChange,
               output_heatmap="MDA231_with_LBH_KDvsMDA231_with_231_siCtrl",
               out_dir= Out_dir)



############################### merge all raw (with statistics and annotation) and normalized reads as well as annotation 
head(PBS_vs_IL2)
dim(PBS_vs_IL2)

head(normSF.Karoline_Briegel.wAnno)
dim(normSF.Karoline_Briegel.wAnno)

normSF.Karoline_Briegel_df  <- normSF.Karoline_Briegel.wAnno[, c("ensembl_gene_id","IL2_1","IL2_2","IL2_3", "PBS_1","PBS_2","PBS_3")] 
head(normSF.Karoline_Briegel_df)

stat_annoted_raw_and_normalized_reads <- merge(PBS_vs_IL2, normSF.Karoline_Briegel_df, by= "ensembl_gene_id")   # "row.names" or the number 0 specifies the row names
stat_annoted_raw_and_normalized_reads <- stat_annoted_raw_and_normalized_reads[order(stat_annoted_raw_and_normalized_reads$padj), ]

head(stat_annoted_raw_and_normalized_reads)
dim(stat_annoted_raw_and_normalized_reads)

write.table(stat_annoted_raw_and_normalized_reads, file.path(Out_dir, "stat_annoted_raw_and_normalized_reads.xls"), row.names=F, sep = "\t")
###############################


```
