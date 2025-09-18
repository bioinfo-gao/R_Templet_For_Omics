require(ISLR)
names(Smarket)
summary(Smarket)
?(Smarket)
pairs(Smarket, col=Smarket$Direction)

glm.fit=glm(data = Smarket, Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, family = binomial)
summary(glm.fit)
