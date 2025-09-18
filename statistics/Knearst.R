require(ISLR)
require(class)
?knn
ls()
objects(2) # from the above ls()
attach(Smarket)
names(Smarket)
summary(Smarket)
head(Smarket)
?(Smarket)
#pairs(Smarket, col=Smarket$Direction)

## 
Xlag=cbind(Lag1, Lag2)
train =Year<2005
Xlag[train,]
Xlag[!train,]
head(Direction)

head(Smarket[, Direction])
DDirection=Smarket[train, c("Year","Direction")] # DDirection= Smarket[, c(Year, Direction), drop = F]

head(DDirection)


dim(Xlag[train,])
DDDirection = DDirection[, "Direction"]

DDDirection=Smarket[, c("Year","Direction")] # DDirection= Smarket[, c(Year, Direction), drop = F]

#knn.pred=knn(Xlag[train, ], Xlag[!train,], DDirection[ Direction] , k=1)
knn.pred=knn(Xlag[train, ], Xlag[!train,], DDDirection , k=1)

knn.fit
summary(knn.fit)
plot(knn.fit)
Smarket.2005= subset (Smarket, Year==2005)
knn.pred=predict(knn.fit,  Smarket.2005)
knn.pred[,]
dim(knn.pred)
summary(knn.pred)
class(knn.pred)
data.frame(knn.pred)[1:5,]
