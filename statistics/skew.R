install.packages("moments")
require(moments)
require(ggplot2)

sourceAddress <- "http://www.olympic.org/olympic-results/london-2012/athletics/long-jump-m"

longjump.df <- read.csv(header=TRUE, sep=",", text="
rank,name,country,distance
1,Mauro Vinicius DA SILVA,BRA,8.11 
2,Marquise GOODWIN,USA,8.11
3,Aleksandr MENKOV,RUS,8.09
4,Greg RUTHERFORD,GBR,8.08
5,Christopher TOMLINSON,GBR,8.06
6,Michel TORNEUS,SWE,8.03
7,Godfrey Khotso MOKOENA,RSA,8.02
8,Will CLAYE,USA,7.99
9,Mitchell WATT,AUS,7.99,
10,Tyrone SMITH,BER,7.97,
11,Henry FRAYNE,AUS,7.95,
12,Sebastian BAYER,GER,7.92,
13,Christian REIF,GER,7.92,
14,Eusebio CACERES,ESP,7.92,
15,Aleksandr PETROV,RUS,7.89,
16,Sergey MORGUNOV,RUS,7.87,
17,Mohammad ARZANDEH,IRI,7.84,
18,Ignisious GAISAH,GHA,7.79,
19,Damar FORBES,JAM,7.79,
20,Jinzhe LI,CHN,7.77,
21,Raymond HIGGS,BAH,7.76,
22,Alyn CAMARA,GER,7.72,
23,Salim SDIRI,FRA,7.71,
24,Ndiss Kaba BADJI,SEN,7.66,
25,Arsen SARGSYAN,ARM,7.62,
26,Povilas MYKOLAITIS,LTU,7.61,
27,Stanley GBAGBEKE,NGR,7.59,
28,Marcos CHUVA,POR,7.55,
29,Louis TSATOUMAS,GRE,7.53,
30,Stepan WAGNER,CZE,7.50,
31,Viktor KUZNYETSOV,UKR,7.50,
32,Luis RIVERA,MEX,7.42,
33,Ching-Hsuan LIN,TPE,7.38,
33,Supanara SUKHASVASTI N A,THA,7.38,
35,Boleslav SKHIRTLADZE,GEO,7.26,
36,Xiaoyi ZHANG,CHN,7.25,
37,Mohamed Fathalla DIFALLAH,EGY,7.08,
38,Roman NOVOTNY,CZE,6.96,
39,George KITCHENS,USA,6.84,
40,Vardan PAHLEVANYAN,ARM,6.55,
NA,Luis MELIZ,ESP,NA,
NA,Irving SALADINO,PAN,NA")

roundedSkew <- signif(skewness(longjump.df$distance, na.rm=TRUE), 3)

longjump.df.back = longjump.df
longjump.df$distance = (longjump.df.back$distance)^2
longjump.df

ggplot(longjump.df, aes(x=distance)) + 
    xlab("Distance in metres") +
    ggtitle("London 2012 Men's Long Jump qualifying round results") +
    geom_rug(size=0.8) + 
    geom_density(fill="steelblue") +
    annotate("text", x=7.375, y=0.0625, colour="white", label=paste("Source:", sourceAddress), size=3) +
    #annotate("rect", xmin = 6.25, xmax = 7.25, ymin = min(), ymax = 1.125, fill="white") +
    #annotate("rect", xmin = 6.25, xmax = 7.25, fill="white") +
    scale_y_continuous(expand = expansion(0,0.1)) +
    annotate("text", x=6.75, y=1, colour="black", label="Best jump in up to 3 attempts") +
    annotate("text", x=6.75, y=.875, colour="black", label="42 athletes competed") +
    annotate("text", x=6.75, y=.75, colour="black", label="2 athletes had no legal jump") +
    annotate("text", x=6.75, y=.625, colour="black", label=paste("Skewness = ", roundedSkew))

ggplot(longjump.df, aes(x=distance)) + 
    xlab("Distance in metres") +
    ggtitle("London 2012 Men's Long Jump qualifying round results") +
    geom_rug(size=0.8) + 
    geom_density(fill="steelblue") +
    annotate("text", x=7.375, y=0.0625, colour="white", label=paste("Source:", sourceAddress), size=3) +
    annotate("rect", xmin = 6.25, xmax = 7.25, ymin = 0.5, ymax = 1.125, fill="white") +
    annotate("text", x=6.75, y=1, colour="black", label="Best jump in up to 3 attempts") +
    annotate("text", x=6.75, y=.875, colour="black", label="42 athletes competed") +
    annotate("text", x=6.75, y=.75, colour="black", label="2 athletes had no legal jump") +
    annotate("text", x=6.75, y=.625, colour="black", label=paste("Skewness = ", roundedSkew))


# Results of the top twelve who qualified for the Final were closer to symmetric
skewness(longjump.df$distance[1:12])
# -0.1248782

# Results in the Final (some had 3 jumps, others 6) were only slightly negatively skewed
skewness(c(8.31, 8.16, 8.12, 8.11, 8.10, 8.07, 8.01, 7.93, 7.85, 7.80, 7.78, 7.70))
# -0.08578357

# Compare to Seoul 1988 Heptathlon
install.packages("HSAUR")
require(HSAUR)
skewness(heptathlon)
