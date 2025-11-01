require(ISLR)
require(MASS)
names(Smarket)
summary(Smarket)
?(Smarket)

#pairs(Smarket, col=Smarket$Direction)

## Linear Discrimiant Analysis  
lda.fit=lda(data = Smarket, Direction ~ Lag1 + Lag2, subset =Year<2005)
lda.fit
summary(lda.fit)
plot(lda.fit)
Smarket.2005= subset (Smarket, Year==2005)
lda.pred=predict(lda.fit,  Smarket.2005)
lda.pred[,]
dim(lda.pred)
summary(lda.pred)
class(lda.pred)
data.frame(lda.pred)[1:5,]
