library(corrplot)
library(ggplot2)
library(MASS)
library(leaps)
library(glmnet) 
library(caret)
library(dplyr)
library(readr)
library(nnet)
library(car)
library(base)
library(mcmcplots)


wine_data <-read.csv("winequality-red.csv", TRUE, ";")
attach(wine_data)
cor(quality, density)
model <- lm(quality~density)
summary(model)

str(wine_data)

## histogram
par(mfrow=c(3,4))
hist(fixed.acidity, col="indianred3")
hist(volatile.acidity, col="indianred3")
hist(citric.acid, col="indianred3")
hist(residual.sugar, col="indianred3")
hist(chlorides, col="indianred3")
hist(free.sulfur.dioxide, col="indianred3")
hist(total.sulfur.dioxide, col="indianred3")
hist(density, col="indianred3")
hist(pH, col="indianred3")
hist(sulphates, col="indianred3")
hist(alcohol, col="indianred3")
barplot(table(quality), col=c("indianred4", "indianred3", "indianred2", "indianred1", "indianred", "indianred4"))
mtext("Quality", side=1, outer=F, line=2, cex=0.5)

## correlation plot
corrplot(cor(wine_data))

## scatterplot
par(mfrow=c(2,4))
plot(fixed.acidity, citric.acid, main="Scatterplot", xlab="fixed.acidity", ylab="citric.acid", 
     ylim = c(0, 1.2), xlim = c(4, 16), cex=0.5, col=8)
abline(lm(citric.acid~fixed.acidity), col=2)

plot(fixed.acidity, density, main="Scatterplot", xlab="fixed.acidity", ylab="density", 
     xlim = c(4, 16), cex=0.5, col=8)
abline(lm(density~fixed.acidity), col=2)

plot(fixed.acidity, pH, main="Scatterplot", xlab="fixed.acidity", ylab="pH", col=8)
abline(lm(pH~fixed.acidity), col=2)

plot(free.sulfur.dioxide, total.sulfur.dioxide, main="Scatterplot", xlab="free.sulfur.dioxide", ylab="total.sulfur.dioxide", col=8)
abline(lm(total.sulfur.dioxide~free.sulfur.dioxide), col=2)

plot(citric.acid, volatile.acidity, main="Scatterplot", xlab="citric.acid", ylab="volatile.acidity", cex=0.6, col=8)
abline(lm(volatile.acidity~citric.acid), col=2)

plot(citric.acid, pH, main="Scatterplot", xlab="citric.acid", ylab="pH", cex=0.6, col=8)
abline(lm(pH~citric.acid), col=2)

plot(alcohol, density, main="Scatterplot", xlab="alcohol", ylab="density", cex=0.6, col=8)
abline(lm(density~alcohol), col=2)

plot(pH, quality, main="Scatterplot", xlab="pH", ylab="quality", cex=0.6, col=8)
abline(lm(quality~pH), col=2)

plot(alcohol, quality, main="Scatterplot", xlab="alcohol", ylab="quality", cex=0.6, col=8, type="l")
abline(lm(quality~alcohol), col=2)

## train data set and test data set
train.set <- seq(1,1120)  ## 80%
train_data <- wine_data[train.set,]

test.set <- seq(1121,1599)  ## 30%
test_data <- wine_data[test.set,]

model <- lm(quality~scale(volatile.acidity)+scale(citric.acid)+scale(residual.sugar)+scale(chlorides)+
              scale(free.sulfur.dioxide)+scale(total.sulfur.dioxide)+scale(pH)+scale(sulphates)+scale(alcohol), data=train_data)

vif(model)
summary(model)
AIC (model)
model1 <- lm(quality~scale(volatile.acidity)+scale(citric.acid)+scale(chlorides)+
              scale(total.sulfur.dioxide)+scale(sulphates)+scale(pH)+
               scale(alcohol), data = train_data)

summary(model1)
AIC(model1)
BIC(model1)
vif(model1)
## standardized variables
x1 <- scale(train_data$volatile.acidity)
x2 <- scale(train_data$citric.acid)
x3 <- scale(train_data$residual.sugar)
x4 <- scale(train_data$chlorides)
x5 <- scale(train_data$free.sulfur.dioxide)
x6 <- scale(train_data$total.sulfur.dioxide)
x7 <- scale(train_data$sulphates)
x8 <- scale(train_data$pH)
x9 <- scale(train_data$alcohol)

x1_t <- scale(test_data$volatile.acidity)
x4_t <- scale(test_data$chlorides)
x6_t <- scale(test_data$total.sulfur.dioxide)
x7_t <- scale(test_data$sulphates)
x9_t <- scale(test_data$alcohol)

X_test <- cbind(1, x1_t, x4_t, x6_t, x7_t, x9_t)
Y_obs <- test_data$quality

X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9)


y <- as.matrix(train_data$quality)


lpy.X <- function (y, X, g=length(y),
  nu0=1, s20=try(summary(lm(y~X-1,data=train_data))$sigma^2, silent=TRUE)){
  n<-dim(X)[1] ; p<-dim(X)[2]
  if (p==0){s20 <- mean(y^2)}
  SSRg <- ifelse(p>0, 
                 t(y)%*%(diag (1,nrow=n)-(g/( g+1)) * X%*%solve(t(X)%*%X)%*%t(X))%*%y,
                 t(y)%*%(diag (1,nrow=n)-0)%*%y)
  -.5*(n*log(pi)+p*log(1+g)+(nu0+n)*log(nu0*s20+SSRg)-nu0*log(nu0*s20))+lgamma((nu0+n)/2)-lgamma(nu0/2)
}

z <- rep(1, dim(X)[2])
lpy.c <- lpy.X(y,X[, z==1, drop=FALSE])
S <- 10000
Z <- matrix(NA,S,dim(X)[2])

## Gibbs sampling 
for (s in 1:S){
  for (j in sample(1:dim(X)[2])){
    zp<-z ; zp[j]<-1-zp[j]
    lpy.p<- lpy.X(y,X[,zp==1, drop=FALSE])
    r<- (lpy.p-lpy.c)*(-1)^(zp[j]==0)
    z[j]<-rbinom(1,1,1/(1+exp(-r)))
    if(z[j]==zp[j]) {lpy.c<-lpy.p}
  }
  Z[s,]<-z
  Xnew <- X[, zp==1, drop=FALSE]
  
  n<-length(y)
  p<-dim(Xnew)[2]
  g <- length(y)
  
  nu0=1
  s20=try(summary(lm(y~Xnew-1,data=train_data))$sigma^2, silent=TRUE)
  
  nu <- (nu0+n)/2
  Hg <-  (g/( g+1)) * Xnew%*%solve(t(Xnew)%*%Xnew)%*%t(Xnew) 
  tau <- (nu0*s20 + t(y)%*%(diag (1,nrow=n)-Hg)%*%y)/2
  
  s2 <- 1/rgamma(S, nu, tau)
  
  beta.ols<- solve(t(Xnew)%*%Xnew)%*%t(Xnew)%*%y
  mu <-  (g/( g+1)) *beta.ols
  Sigma <- (g/( g+1)) * solve(t(Xnew)%*%Xnew) 
  betastd <- matrix(rnorm(S*p, 0, sqrt(s2)), S, p)
  beta <- t(t(betastd%*%chol(Sigma))+c(mu))
    
  print(s)
}

sum(Z[,1])/10000
sum(Z[,2])/10000
sum(Z[,3])/10000
sum(Z[,4])/10000
sum(Z[,5])/10000
sum(Z[,6])/10000
sum(Z[,7])/10000
sum(Z[,8])/10000
sum(Z[,9])/10000
sum(Z[,10])/10000

par(mfrow=c(2,2))
x<- t(as.matrix(c(1,2,3,4,5,6,7,8,9,10)))
Y<- t(as.matrix(c(1,1,0.0551,0.0588, 0.9368, 0.056,1,1, 0.3985,1)))
plot(x,Y, xlab="regressors index", ylab="P(zj=1|y,X)", ylim=c(0,1.2))
abline(h=0.8)
title(main = "Estimated Posterior Probabilities of z=1")

mcmcplot(Xnew)
mcmcplot(beta)
mcmcplot(s2)

beta_new <- c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,5]), median(beta[,6]), median(beta[,7]))
beta_new <- as.matrix(beta_new)
X_new <- Xnew[,c(1:3,5:7)]
Y <- X_new%*%beta_new+rnorm(1120,0,sqrt(median(s2)))
plot(y-Y, ylab="Residuals",  ylim=c(-10,10))
abline(h=0, col="red")
qqnorm(y-Y)
qqline(y-Y)
plot(Y)
points(y,col="red") 



### prediction

Y_predict <- X_test%*%beta_new+rnorm(479,0,sqrt(median(s2)))

## residual plot
plot(Y_obs-Y_predict, ylab="Residuals", ylim=c(-10,10))
abline(h=0, col="red")

## average squared predictive error
mean((Y_obs-Y_predict)^2)

qqnorm(Y_obs-Y_predict)
qqline(Y_obs-Y_predict)

plot(Y_predict, Y_obs-Y_predict, ylab="Residuals", xlab="Predicted Values")
plot(Y_predict, Y_obs)
abline(a=0,b=1, col="red")

## quantile
quantile(beta[,1],c(0.025,0.975))
quantile(beta[,2],c(0.025,0.975))
quantile(beta[,3],c(0.025,0.975))
quantile(beta[,5],c(0.025,0.975))
quantile(beta[,6],c(0.025,0.975))
quantile(beta[,7],c(0.025,0.975))