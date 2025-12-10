x <- seq(10)
x
?read.table
?read.table
?read.table
?read.csv2
?read.csv2
?UCBAdmissions
str(UCBAdmissions)
UCBAdmissions
dim(UCBAdmissions)
plot(UCBAdmissions)
margin.table(UCBAdmissions,1)
margin.table(UCBAdmissions,2)
margin.table(UCBAdmissions,3)
admit.dept <- margin.table(UCBAdmissions,3)
admit.dept
barplot(admit.dept)
prop.table(admit.dept)
round(prop.table(admit.dept),2)
admit1 = as.data.frame.table(admit.dept)
admit1 
admit2 = as.data.frame.table(UCBAdmissions)
admit2 
admit2[1:5,3:4] 

# ============

x = c(12,4,45,5,68,45,58)
#browseURL("http://research.stowers−institute.org/efg/R/Color/Chart")
browseURL("https://www.scribd.com/document/410908677/Chart-of-R-Colors")
colors()

col2rgb("navyblue")
barplot(x, col= "navyblue")
barplot(x, col= "slategray3")
barplot(x, col= colors()[602])
barplot(x, col= colors()[102])
barplot(x, col= rgb(159,182,205,max=255))
barplot(x, col= rgb(0.54,0,0))
barplot(x, col= "#9FB6CF")
barplot(x, col= "#9FB6CF")

help(package  = colorspaces)
palette()
barplot(x, col= 1:6)
barplot(x, col= rainbow(6))
barplot(x, col= heat.colors(6))
barplot(x, col= terrain.colors(6))
barplot(x, col= topo.colors(6))
barplot(x, col= cm.colors(6))
