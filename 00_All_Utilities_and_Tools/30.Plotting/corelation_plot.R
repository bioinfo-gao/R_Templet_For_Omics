#devtools::install_github("Hy4m/linkET")
# install.packages("vegan")

library(ggplot2)
library(linkET) # 核心包。用于相关性和热图
library(vegan)  # Mantel检验和距离计算
library(dplyr)  
library(RColorBrewer)#」色漸变自定义)

#mtcars的相关(Pearson方法,默认)
cor_mtcars <-  correlate(mtcars)

# ?as_md_tbl # linkET
#格式,后续给图
as_md_tbl(cor_mtcars)
data(varespec)
data(varechem)
head(varespec)
head(varechem)

                    #测值:弱/中/强相关 
                    # 阀值:高度/中度/无显著 
mantel_results <- mantel_test(varespec, 
                              varechem,
                              spec_select = list(Group1 =  1:10,
                                                 Group2 = 11:25,
                                                 Group3 = 26:35,
                                                 Group4 = 36:44)) %>% 
    mutate(rd = cut(r, breaks = c(-Inf,   0.3,  0.5, Inf), labels = c("< 0.3",     "0.3-0.5",  ">= 0.5")),
           pd = cut(p, breaks = c(-Inf, 0.001, 0.01, Inf), labels = c("< 0.001", "0.001-0.01", ">= 0.01"))
           )


head(mantel_results)
(mantel_results)

#计算:使用bray 距离 for spec, euclidean for env

#绘制热图+ Mante线条叠加
qcorrplot(correlate(varechem), type="lower", diag = FALSE) +  # 下三角,无对角线 
    geom_tile() + 
    geom_couple(aes(colour= pd, size = rd), 
                data = mantel_results, 
                curvature =  nice_curvature()) +  # 优雅曲线连接 
    scale_fill_gradientn(colours = brewer.pal(11, "BrBG")) + #BrBG方案(薪) 
    scale_size_manual(values = c(0.8, 1.5,3)) +           # 更细腻 
    scale_colour_manual(values = c("#004488",  "#DDAA33" , "#BB5566")) + #自定义配色,深师全黄红
    guides(size = guide_legend(title="Mantel's r", override.aes = list(colour = "grey50"), order = 2),
           colour = guide_legend(title="Mantel's p", override.aes = list(size = 4), order = 1),
            fill = guide_colorbar(title="Pearson's r", order = 3))


## <<<<<<<<<<<<====================== 注意：R　和P　事实上只有两个区间，所有只能作图如此














