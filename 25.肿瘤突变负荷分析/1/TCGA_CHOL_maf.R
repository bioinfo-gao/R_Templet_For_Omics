
#这个函数有两个参数，
#@metadata是从TCGA数据下载的sample sheet
#@path是保存maf文件的路径

merge_maf <- function(metadata, path){
  #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)

  message ('############### Merging maf data ################\n',
           '### This step may take a few minutes ###\n')
    #通过lapply循环去读每一个样本的maf，然后通过rbind合并成矩阵，按行来合并
    #colClasses指定所有列为字符串
    mafMatrix <- do.call("rbind", lapply(filenames, function(fl) 
      read.table(gzfile(fl),header=T,sep="\t",quote="",fill=T,colClasses="character")))
    return (mafMatrix)
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



#读入maf的sample sheet文件
metaMatrix.maf=read.table("maf_sample_sheet.tsv",sep="\t",header=T)
#替换.为下划线，转换成小写，sample_id替换成sample
names(metaMatrix.maf)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.maf))))
#删掉最后一列sample_type中的空格
metaMatrix.maf$sample_type=gsub(" ","",metaMatrix.maf$sample_type)

#删掉重复的样本
metaMatrix.maf <- FilterDuplicate(metaMatrix.maf)


#调用merge_maf函数合并maf的矩阵
maf_value=merge_maf(metadata=metaMatrix.maf, 
                     path="maf_data"
                     )
#查看前三行前十列
maf_value[1:3,1:10]


#保存合并后的maf文件
write.table(file="combined_maf_value.txt",maf_value,row.names=F,quote=F,sep="\t")


