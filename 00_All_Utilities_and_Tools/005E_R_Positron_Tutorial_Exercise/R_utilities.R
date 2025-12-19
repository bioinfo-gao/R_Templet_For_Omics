plot(1:26, 1:26, pch=1:26, cex=2)

## Load the data and fit a multiple regression model
ps <- read.csv("https://remiller1450.github.io/data/Salaries.csv")
m <- lm(yrs.service ~ yrs.since.phd + rank, data = ps)

## Hypothesis tests
summary(m)

#install.packages("ISLR")
library("MASS")
library("ISLR")
?Boston
?Carseats

attach("ISLR")
data(Carseats) # 2024 must load for fix () fix(Carseats) ,otherwise it will fix a new funtion
fix(Carseats)
summary(Carseats)

