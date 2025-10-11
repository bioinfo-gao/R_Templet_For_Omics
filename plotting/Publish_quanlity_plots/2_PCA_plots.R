# Load the USArrests dataset
data(USArrests)

# View the first few rows of the dataset
head(USArrests)

# Standardize the data (important for PCA, especially when variables have different scales)
scaled_data <- scale(USArrests)

# Perform PCA using prcomp()
# The 'scale = TRUE' argument in prcomp() will automatically scale the data,
# so explicitly scaling beforehand is optional if using this argument.
# However, scaling manually can offer more control or be useful for other analyses.
pca_result <- prcomp(USArrests, scale = TRUE)

# View the PCA results
print(pca_result)

# Get a summary of the PCA results, including proportion of variance explained
summary(pca_result)

# Extract the principal component scores (the "PCA dataset")
pca_scores <- as.data.frame(pca_result$x)

# View the first few rows of the PCA scores
head(pca_scores)

# You can also access the loadings (eigenvectors)
pca_loadings <- pca_result$rotation
print(pca_loadings)



#载入ggplot2绘图包；
library(ggplot2)
#安装ggh4x包；
#install.packages("ggh4x")
library(ggh4x)


dt = (pca_scores)


dt$class =  c("East", "West", "West", "West", "West", "West", "East", "East", "East", "East",
              "West", "West", "East", "East", "West", "West", "East", "East", "East", "East",                "East", "West", "West", "West", "West", "West", "East", "East", "East", "West",
              "West", "West", "East", "East", "East", "East", "East", "West", "West", "West", 
              "West", "West", "East", "East", "East", "East", "East", "West", "West", "West")

table(dt$class) 
dt[ , "class", drop = FALSE ] 
 
#绘制散点图；
?ggplot2
#p1 <- ggplot(dt,aes(x=PC1,y=PC2, color= rownames(dt)))+ 
p1 <- ggplot(dt,aes(x=PC1,y=PC2, color= class )) + 
  geom_point()

#预览初始绘图效果；
p1

#自定义分组颜色、透明度（对应色号的后两位）、图例标题；
p2 <- p1+scale_color_manual(name="",
                            values = c("#e26fffbb","#8e9cffbb"))
#预览绘图效果；
p2

#自定义坐标轴范围和刻度；
#通过guide = "axis_minor"添加小刻度；
p3 <- p2+scale_x_continuous(position = "bottom",
                            guide = "axis_minor",
                            expand=expansion(add = c(0.7,0)),
                            limits=c(-10,5),
                            breaks = c(-10,-5,0, 5),
                            label = c("-10","-5", "0", "5"))+
  scale_y_continuous(position = "left",
                     guide = "axis_minor",
                     expand=expansion(add = c(0.5,0)),
                     limits=c(-7.5,5),
                     breaks = c(-7.5,-5,-2.5,0,2.5,5),
                     label = c("-7.5","-5", "-2.5", "0","2.5","5"))

p3




#自定义主题；
#包括添加小刻度、自定义刻度长短等；
mytheme <- theme(panel.grid = element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 legend.position = "top",
                 axis.line = element_line(colour = "grey30"),
                 axis.ticks.length = unit(1.8, "mm"),
                 ggh4x.axis.ticks.length.minor = rel(0.6))
p4 <- p3+mytheme
p4




#调整坐标轴的样式，调整坐标轴线的长度；
#比较遗憾的是当前版本无法保留小刻度,与guides(x = "axis_minor", y = "axis_minor")冲突；
#自定义图例分组图形（点）的大小；
p5 <- p4+guides(x = guide_axis_truncated(trunc_lower = -10,
                                         trunc_upper = 5),
                y = guide_axis_truncated(trunc_lower = -7.5,
                                         trunc_upper = 5),
                color = guide_legend(override.aes = list(size=4)))

p5



#综上技巧，尝试添加质心线；
p6 <- ggplot(dt,aes(x=PC1,y=PC2,color=class))+
  stat_centroid(aes(xend = PC1, yend = PC2, colour = class),
                geom = "segment", crop_other = F,
                alpha=0.3,size = 1,show.legend = F)+
  geom_point(size=2,stroke=1,alpha=1,
             fill="white",shape = 21,show.legend = T)+
  scale_color_manual(name="",
                     values = c("#e26fff","#8e9cff"))+
  scale_x_continuous(expand=expansion(add = c(0.7,0.7)),
                     limits=c(-10,5))+
  scale_y_continuous(expand=expansion(add = c(0.5,0.5)),
                     limits=c(-7.5,5))+
  guides(x = "axis_truncated",y = "axis_truncated")+mytheme

p6




#绘制实心散点图；
p7 <- ggplot(dt,aes(x=PC1,y=PC2,fill=class))+
  stat_centroid(aes(xend = PC1, yend = PC2, colour = class),
                geom = "segment", crop_other = F,
                alpha=0.3,size = 1,show.legend = F)+
  geom_point(size=3,alpha=0.7,
             color="white",shape = 21,show.legend = T)+
  scale_color_manual(name="",
                     values = c("#FF9999","#c77cff"))+
  scale_fill_manual(name="",
                    values = c("#FF9999","#c77cff"))+
  scale_x_continuous(expand=expansion(add = c(0.7,0.7)),
                     limits=c(-10,5))+
  scale_y_continuous(expand=expansion(add = c(0.5,0.5)),
                     limits=c(-7.5,5))+
  guides(x = "axis_truncated",y = "axis_truncated")+mytheme

p7





  
  
  
  
