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