myda<-data.frame( 
  disease=c("流感","感染性腹泻","手足口病","新冠肺炎","病毒性肝炎","急性出血\n性结膜炎",  
            "梅毒","肺结核","淋病","流行性\n腮腺炎","百日咳","艾滋病","登革热","水痘","猩红热",  
            "伤寒","布病","痢疾"),
  ynn=c(2000,1500,1000,920,800,720,630,510,420,360,290,210, 
        150,88,68,47,30,9))

library(viridis)
ggplot(data = myda,aes(x=reorder(disease,-ynn),y=ynn,fill=disease))+ 
  geom_bar(width = 0.8,stat = "identity")+  coord_polar(theta = "x",start=0)+  
  ylim(-200,2000)+scale_fill_viridis(option="D",discrete=T)


ggplot(data = myda,aes(x=reorder(disease,-ynn),y=ynn,fill=disease))+  
  geom_bar(width = 0.8,stat = "identity")+  coord_polar(theta = "x",start=0)+  
  ylim(-200,2000)+  scale_fill_viridis(option = "D",discrete = T)+  
  theme_minimal()+xlab(" ")+ylab(" ")+ 
  labs(    title = "欢迎关注公众号:R语言与医学生", 
           subtitle = paste(  "视频教学可见B站:R语言与医学生。",      sep = "\n"    ),   
           caption = "2024.01.06")+theme(legend.position="none")

# 2.7添加文本
# 我们为图形添加发病数，最简单的就是通过geom_text函数实现。

ggplot(data = myda,aes(x=reorder(disease,-ynn),y=ynn,fill=disease))+  
  geom_bar(width = 0.8,stat = "identity")+  coord_polar(theta = "x",start=0)+  
  ylim(-400,2000)+  scale_fill_viridis(option = "D",discrete = T)+  
  theme_minimal()+xlab(" ")+ylab(" ")+  labs(  
    title = "欢迎关注公众号:R语言与医学生",
    
    subtitle = paste(      "\n公众号:R语言与医学生,致力于R语言代码分享!",  
                           
                           "视频教学可见B站:R语言与医学生。",      sep = "\n"    ),  
    caption = "2024.01.06")+  theme(legend.position = "none")+
  geom_text(aes(label=ynn),size=4,vjust=-0.5,hjust=0.5,color="#8B8878")


####构造标签
label_data<-myda
#计算文本的角度，首先计算有多少个文本
#number_of_bar = 18
library(data.table)
setDT(label_data)
#构造文本
label_data[,new_label:=paste0(disease,ynn,"例")]
label_data$id=1:nrow(label_data)

number_of_bar <- nrow(label_data)
label_data[ , angle:=90 - 360 * (label_data$id-0.5) /number_of_bar]

label_data[,":="(hjust=ifelse(angle<90,1,0),angle1=ifelse(angle<90,angle+180,angle))]

label_data

# Sys.setlocale("LC_ALL", "zh_CN.UTF-8")

ggplot(data = myda,aes(x=reorder(disease,-ynn),y=ynn,fill=disease))+ 
  geom_bar(width = 0.8,stat = "identity")+  coord_polar(theta = "x",start=0)+  
  ylim(-400,2000)+  scale_fill_viridis(option = "D",discrete = T)+ 
  theme_minimal()+xlab(" ")+ylab(" ")+ 
  theme(  
    legend.position = "noine",   
    text = element_text(color = "gray12", family = "Bell MT"),  
    axis.text = element_blank(),    axis.title = element_blank(),  
    panel.grid = element_blank(),    plot.margin = unit(rep(-1,4), "cm")   )+ 
  geom_text(data=label_data, aes(x=id, y=ynn, label=new_label, hjust=hjust),  
            color="black", fontface="bold",  
            alpha=0.6, size=3.5, angle=label_data$angle1,inherit.aes=FALSE)

