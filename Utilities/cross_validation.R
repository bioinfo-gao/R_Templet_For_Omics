require(ISLR)
require(class)
require(boot)
?cv.glm


summary(Auto)

plot(mpg ~ horsepower, data=Auto)

#LOOCV
glm.fit = glm(mpg ~ horsepower, data=Auto)
cv.glm(Auto, glm.fit)$delta

loocv=function(fit){
    h=lm.influence(fit)$h
    mean(residuals(fit)/(1-h)^2)
    
}

loocv(glm.fit)

cv.error=rep(0,5)

degree=1:5

for (d in degree){
    glm.fit=glm( mpg~poly(horsepower, d) , data = Auto)
    cv.error[d]=loocv(glm.fit)
    
}


cv.error
plot(degree, cv.error, type="b")
