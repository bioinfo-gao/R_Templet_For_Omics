# "koenrich.xls" 
# The term "KO enrichment" refers to KEGG Orthology enrichment, a common bioinformatics analysis
# 
# A koenrich.xls file is a standard output from bioinformatics tools like DAVID or those based on the clusterProfiler
#  representing the results of a KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment analysis. 
#  The spreadsheet show which KEGG pathways are statistically over-represented (or "enriched") in a given gene list. 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")  
# browseVignettes("clusterProfiler")

getwd()
#安装所需R包，安装过的无需安装；
#install.packages("ggplot2")
#读入数据；
dt1<-read.table("koenrich.xls",sep = ",",header = T)
head(dt1,16)
head(dt1)
dim(dt1)
colnames(dt1)

#指定纵轴标签顺序,按照输入文件的顺序排序，否则默认按照首字母顺序,同时逆序绘制，保持与表格顺序一致；
dt1$KEGG_A_Class<-factor(dt1$KEGG_A_Class,
                         levels = rev(unique(dt1$KEGG_A_Class)),
                         
                         ordered = TRUE)

dt1$KEGG_B_Class<-factor(dt1$subcategory,
                         levels = rev(unique(dt1$subcategory)),
                         ordered = TRUE)           


# 加载ggplot2包；
library(ggplot2)
#建立数据(Genes.Number, KEGG_B_Class)与图形（点）的映射关系，即确定点的(x,y)坐标，绘制散点图；
p1<-ggplot(dt1, aes(x=Count, #p1<-ggplot(dt1, aes(x=Genes_Number, 
                    y=KEGG_B_Class,
                    fill=KEGG_A_Class,na.rm = FALSE))+
  geom_bar(stat="identity",na.rm = FALSE)
p1

#添加标签,nudge_y，y方向偏移1.5，；
p2<-p1+geom_text(aes(x=Count,
                     y=KEGG_B_Class,
                     label= Count
                     #label=label
                     ),
                 size=2.5,
                 hjust="left",
                 nudge_x= 5)
                 #nudge_x=0.1)
p2

#设置x轴范围，避免文字溢出绘图区；
p3<-p2+scale_x_continuous(limits = c(0, 60),
                          expand=expansion(mult = c(0, .1)))
p3

#设置图例、坐标轴、图表的标题；
p4<-p3+labs(x="Number of Gene",
            y="",
            title="KEGG pathway anotation")
p4

#自定义颜色；
p5<-p4+scale_fill_brewer(palette="Dark2")
p5

#还是默认的配色好看一点，这里直接回看看p4的效果；
#p5 <- p4
#获取颜色；
g <- ggplot_build(p5)
mycol<-g$data[[1]]["fill"]
col<-rev(mycol[,1])
#将A Class对应的颜色设为深黑色；
num <- rev(dt1$Genes_Number)
index <- which(num==0)
col[index] <- "grey10"

#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
mytheme1<-theme(plot.title = element_text(size = rel(1),hjust = 0.5,face = "bold"),
                axis.title = element_text(size = rel(1)),
                axis.text.y = element_text(size=rel(0.85),
                                           colour =col,face = "bold"),
                legend.position = "none",
                plot.margin=unit(x=c(top.mar,right.mar,
                                     bottom.mar,left.mar),
                                 units="inches"))
#查看绘图效果；
p6<-p5+mytheme1
p6




#再换一种主题；
mytheme2<-theme_bw()+theme(plot.title = element_text(size = rel(1),
                                                     hjust = 0.5,face = "bold"),
                           axis.title = element_text(size = rel(1)),
                           axis.text.y = element_text(size=rel(0.85),
                                                      colour =col,face = "bold"),
                           legend.position = "none",
                           panel.grid = element_blank(),
                           plot.margin=unit(x=c(top.mar,right.mar,
                                                bottom.mar,left.mar),
                                            units="inches"))

p7<-p5+mytheme2
p7



# 
# 注意，ggplot2提示“Vectorized input to `element_text()` is not officially supported”，虽然官方不支持，但依然可以用。
# 
# 
# 
# 好啦，今天就分享到这里啦~
  
  
  
