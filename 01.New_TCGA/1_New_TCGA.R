#这个函数有三个参数，
#@metadata是从TCGA数据下载的 ===> sample sheet
#@path是保存样本表达文件的路径

#@data.type是要合并的数据的类型，这里支持RNAseq和miRNAs两种
#@mRNA_expr_type是RNA表达谱的类型，这里支持STARcounts，TPM，FPKM和FPKM-UQ

#@symbol,为T提取基因名字，为F，不提取
#@type,为T提取RNA类型，为F，不提取

setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\01.New_TCGA")#设置工作目录


# FPKM-UQ 是为了解决 FPKM 的文库组成偏差而提出的改进方法，
# 特别在 TCGA (The Cancer Genome Atlas) 等大型项目中被广泛采用。

# 确定上四分位数 (UQ)： 在所有基因的 FPKM 值中，找到第 75 百分位数 (Upper Quartile, UQ) 对应的值。
# 使用 UQ 进行归一化： 最终的 FPKM-UQ 值是基于这个 UQ 值来调整的。
# FPKM-UQ = FPKM 值/ FPKM 值的第 75 百分位数 × 1000


merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    #判断需要合并的是什么数类型，如果是RNAseq执行下面的代码
    if (data.type=='RNAseq') {
        message ('###############    叉叉滴同学   ################\n',
                 '          ### 等我几分钟,这步太吃性能了 ###\n',
                 '          ### 等我几分钟,这步太吃性能了 ###\n',
                 '          ### 等我几分钟,这步太吃性能了 ###\n')
        
        #根据需要的RNAseq表达谱数据类型，提取相应列的数据
        if(      mRNA_expr_type=="STAR"){
            column=4
        }else if(mRNA_expr_type=="TPM"){
            column=7
        }else if(mRNA_expr_type=="FPKM"){
            column=8
        }else if(mRNA_expr_type=="FPKM_UQ"){
            column=9
        }
        
        #通过lapply循环去读每一个样本的表达，然后通过cbind合并成矩阵
        rnaMatrix <- do.call("cbind", 
                             lapply(filenames, 
                                    function(fl) 
                                        read.table(fl,skip=6,sep="\t")[,column]
                             )
        )
        
        #获取第一个文件的第一列作为矩阵的行名
        ensembl     <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
        gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
        type        <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
        
        #排除掉_PAR_Y为后缀的转录本
        index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
        rnaMatrix=rnaMatrix[index,]
        
        #去掉Ensembl ID后面的.和数字，eg.ENSG00000000003.13
        rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
        gene_symbol=gene_symbol[index]
        type=type[index]
        
        #将sample sheet的sample id这一列作为表达矩阵的列名
        colnames(rnaMatrix) <- metadata$sample
        
        
        #统计样本数和基因数
        nSamples = ncol(rnaMatrix)
        nGenes   = nrow(rnaMatrix)
        
        if(RNA_type){
            rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
        }
        
        if(symbol){
            rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
        }
        #输出样本数和基因数
        message (paste('Number of samples: ', nSamples, '\n', sep=''),
                 paste('Number of genes: ', nGenes, '\n', sep=''))
        #返回最后的基因表达矩阵
        return (rnaMatrix)
        
    }else if (data.type=='miRNAs') { #如果需要合并的是miRNA的数据，执行下面代码
        message ('############### Merging miRNAs data ###############\n')
        #利用lapply来读取每个样本miRNA的表达数据，这里需要用到filtermir这个函数
        #filtermir主要提取mature mir的counts数
        mirMatrix <- lapply(filenames, function(fl) filtermir(fl))
        #mirbase是目前人的说有miRNA成熟提的ID号，eg.MIMAT0027618
        mirs <- mirbase$V1
        #cbind合并成矩阵合并成矩阵，这里注意并不是每个样本中都表达所有的miRNA
        #如果不存在某个miRNA,在这一步表达量会用NA表示
        mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                             function(expr) expr[mirs]))
        
        #设置表达矩阵的行名为miRNA的名字，mirbase的第二列就是所有人的miRNA的名字
        rownames(mirMatrix) <- mirbase$V2
        #设置表达矩阵的列名为sample sheet的sample id这一列
        colnames(mirMatrix) <- metadata$sample
        
        #将表达量为NA的地方转换成0
        mirMatrix[is.na(mirMatrix)] <- 0
        
        #统计样本数和miRNA数目
        nSamples = ncol(mirMatrix)
        nGenes = nrow(mirMatrix)
        
        #输出样本数和miRNA的数目
        message (paste('Number of samples: ', nSamples, '\n', sep=''),
                 paste('Number of miRNAs: ', nGenes, '\n', sep=''))
        #返回miRNA的表达矩阵
        return (mirMatrix)
    }else{  #如果data.type不是上面提到的两种，就报错，停止执行
        stop('data type error!')
    }
}

#定义filtermir函数
#@fl参数为所有样本miRNA表达文件的链接
filtermir <- function(fl) {
    #readtable读取文件内容
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    #找到miRNA_region这一列以mature开头的行，eg.mature,MIMAT0000062
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    #将同一个成熟体的所有counts累加起来，作为这个成熟体counts
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    
    #将miRNA_region这一列（mature,MIMAT0000062）按头号分开，取第二列，得到成熟体ID号
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    
    #去掉第一列，只保留miRNA成熟体的counts数
    expr <- expr[,-1]
    #加上对应的成熟体ID作为counts的名字
    names(expr) <- mirs
    #返回过滤之后的miRNA成熟体的counts
    return(expr)
}


#定义去除重复样本的函数FilterDuplicate
FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
        metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
}

#定义过滤样本类型的函数FilterSampleType，只保留PrimaryTumor和SolidTissueNormal这两种样本类型
FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                        c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
        metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
}



#读入RNAseq的sample sheet文件
metaMatrix.RNA=read.table("RNAseq_sample_sheet.tsv",sep="\t",header=T)

#替换.为下划线，转换成小写，sample_id替换成sample

head(metaMatrix.RNA)
dim(metaMatrix.RNA)


names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))

#删掉最后一列sample_type中的空格
metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)

#删掉重复的样本
metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)

#删掉非PrimaryTumor和SolidTissueNormal的样本
#会删掉像Recurrent Tumor这样的样本
metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)

head(metaMatrix.RNA)
dim(metaMatrix.RNA)



#调用merge_TCGA函数合并RNAseq的表达矩阵 << === ~ 10min
RNA_STAR_Counts=merge_TCGA(metadata=metaMatrix.RNA, 
                           path="back/RNAseq", 
                           data.type="RNAseq",
                           mRNA_expr_type="STAR",
                           symbol = T,
                           RNA_type=T
)





### 10 min above << ===
RNA_STAR_Counts[1:10,1:5]
dim(RNA_STAR_Counts)

#保存RNAseq表达矩阵
write.table(file="combined_RNAseq_counts.txt",RNA_STAR_Counts,sep="\t",quote=F)

#TPM
RNA_TPM=merge_TCGA(metadata=metaMatrix.RNA, 
                   path="back/RNAseq", 
                   data.type="RNAseq",
                   mRNA_expr_type="TPM",
                   symbol = T,
                   RNA_type=T
)

RNA_TPM[1:3,1:3] # 2nd 10 min
dim(RNA_TPM)

write.table(file="combined_RNAseq_TPM.txt",RNA_TPM,sep="\t",quote=F)


#FPKM
RNA_FPKM=merge_TCGA(metadata=metaMatrix.RNA, 
                    path="back/RNAseq", 
                    data.type="RNAseq",
                    mRNA_expr_type="FPKM",
                    symbol = T,
                    RNA_type=T
)
RNA_FPKM[1:3,1:3]
write.table(file="combined_RNAseq_FPKM.txt",RNA_FPKM,sep="\t",quote=F)


#FPKM_UQ
RNA_FPKM_UQ=merge_TCGA(metadata=metaMatrix.RNA, 
                       path="back/RNAseq", 
                       data.type="RNAseq",
                       mRNA_expr_type="FPKM_UQ",
                       symbol = T,
                       RNA_type=T
)
RNA_FPKM_UQ[1:3,1:3]
write.table(file="combined_RNAseq_FPKM_UQ.txt",RNA_FPKM_UQ,sep="\t",quote=F)



###############################################
#combine miRNAs
###############################################
#包含人的所有2652个miRNA的名字
load("mirbase.rds")

#读入miRNA seq的sample sheet文件
metaMatrix.MIR=read.table("miRNAs_sample_sheet.tsv",sep="\t",header=T)
#替换.为下划线，转换成小写，sample_id替换成sample
names(metaMatrix.MIR)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.MIR))))
#删掉最后一列sample_type中的空格
metaMatrix.MIR$sample_type=gsub(" ","",metaMatrix.MIR$sample_type)

#删掉重复的样本
metaMatrix.MIR <- FilterDuplicate(metaMatrix.MIR)
#删掉非PrimaryTumor和SolidTissueNormal的样本
#会删掉像Recurrent Tumor这样的样本
metaMatrix.MIR <- FilterSampleType(metaMatrix.MIR)



#调用merge_TCGA函数合并miRNA的表达矩阵
MIR_Counts=merge_TCGA(metadata=metaMatrix.MIR, 
                      path="miRNAs", 
                      data.type="miRNAs"
)
MIR_Counts[1:3,1:3]
#保存RNAseq表达矩阵
write.table(file="combined_miRNAs_counts.txt",MIR_Counts,sep="\t",quote=F)

