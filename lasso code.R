########################################################################################################
#
#   MATH 569 - Final Project
#   Lasso Implementation by - Ankush Agrekar, Divya Vasireddy, Zhiwei Zhang, Jingyuan He
#
########################################################################################################

lassoImp <- function (X, Y, s){
  require ("quadprog")

  p = ncol(X)
  XtX = 2*t(X)%*%X
  dvec = 2*t(X)%*%Y
  Ge = as.matrix(expand.grid(rep(list (c(-1,1)), p)))
  b0 = rep (-s, 2^p)
  bLassoImp = solve.QP (XtX, dvec, -t(Ge), b0)
  return (round(bLassoImp$solution,digits = 6))
}

standardize_data <- function(dataset){
  Xmean <- apply(dataset,2,mean)                # calculate the mean of each column
  temp  <- sweep(dataset,2,Xmean,"-")           # substract the corresponding mean from each element
  Xstdd <- apply(temp,2,var)                    # calculate the variance
  temp  <- sweep(temp,2,sqrt(Xstdd),"/")        # standardized the data to have mean=0 and sd=1
  Data.std  <- as.matrix(temp)
  return(Data.std)
} 

calc_olsbeta <- function(X,Y){
  XtX <- t(X) %*% X
  XtY <- t(X) %*% Y
  olsbeta <- solve(XtX, XtY)
  return(olsbeta)
}

calc_se_Z <- function(X,Y,b){
  yhat <- X %*% b
  errors <- Y - yhat
  MSE <- sum(errors^2)/(nrow(X) - ncol(X))                       # estimate of sigma-squared
  VarCovar <- MSE * chol2inv(chol(t(X) %*% X))                    # variance covariance matrix
  StdErr <- sqrt(diag(VarCovar))
  Zscore <- b / StdErr
  result <- cbind(b, StdErr, Zscore)
  colnames(result) <- c("Coefficients", "Std. Err.", "Z-Score")
  return(result)
}

require(lasso2)
data('Prostate')
Prostate<-Prostate
n <- nrow(Prostate)
p <- ncol(Prostate)

X <- standardize_data(Prostate[,1:8])           # do not pass the response variable
Y <- Prostate[,9]
mX <- cbind(constant = 1, X) 
olsbeta.hat<- calc_olsbeta(mX,Y)

t <- 0.44*sum(abs(olsbeta.hat[-1]))
lasso.calcBeta <- lassoImp(X,as.matrix(Y),t)

mX <- cbind(constant = 1, X)                    # include the column of 1s to include the intercept
beta.new <- rbind(olsbeta.hat[1],matrix(lasso.calcBeta,length(lasso.calcBeta),1))  # including beta0 in the vector of beta

summary.lasso <- calc_se_Z(mX, Y, beta.new)
summary.lasso

########################################################################################################
#
#   comparing with the available function l1ce from lasso2 package
#
########################################################################################################
require(lasso2)
df <- as.data.frame(cbind(X,Y))
names(df)[ncol(df)] <- c("lpsa")
              
lasso.l1ce <- l1ce(lpsa~.,data =df,bound = 0.44)
lasso.l1ce$coefficients

#############################################################################################
#                   LassoImp                                          l1ce
#     Coefficients  Std. Err.   Z-Score                         Coefficients
#[1,]     2.478387 0.07979754 31.058438                           2.47838688
#[2,]     0.558766 0.11496136  4.860468                           0.55876613
#[3,]     0.097000 0.09366990  1.035551                           0.09699974
#[4,]     0.000000 0.09228224  0.000000                           0.00000000
#[5,]     0.000000 0.09407535  0.000000                           0.00000000
#[6,]     0.155587 0.11220755  1.386600                           0.15558749
#[7,]     0.000000 0.14118132  0.000000                           0.00000000
#[8,]     0.000000 0.12614998  0.000000                           0.00000000
#[9,]     0.000000 0.13833647  0.000000                           0.00000000
#

########################################################################################################
#
#   s.e. estimation using BootStrapping
#
########################################################################################################
lasso.calcBeta.boot <- matrix(0,200,ncol(X)+1, byrow = TRUE)
colnames(lasso.calcBeta.boot) <- c("intercept", colnames(X))
lasso.se.boot <- matrix(0,200,ncol(X)+1,byrow = TRUE)

for(iter in 1:200){
  idx <- sample(n, (2*n)/3,replace = TRUE)
  trainData <- X[idx,]
  trainResp <- Y[idx]
  
  mX.boot <- cbind(constant = 1, trainData)
  boot.olsbeta <- calc_olsbeta(mX.boot,trainResp)
  t <- 0.44*sum(abs(olsbeta.hat[-1]))
  lasso.calcBeta <- lassoImp(trainData,as.matrix(trainResp),t)
  
  beta.new <- rbind(round(boot.olsbeta[1],digits = 6),matrix(lasso.calcBeta,length(lasso.calcBeta),1))
  lasso.calcBeta.boot[iter,] <- beta.new
  lasso.se.boot[iter,] <- calc_se_Z(mX.boot, trainResp, beta.new)[,2]
}

mean.se.boot <- apply(lasso.se.boot,2,mean)

par(mfcol=c(1,p))
boxplot(lasso.calcBeta.boot[,2:ncol(X)+1])



########################################################################################################
#
#   Simulation part with mtcars.
#
########################################################################################################
library(leaps)
library(MASS)

carsData <- mtcars
ncar <- nrow(carsData)
pcar <- ncol(carsData)

Xcars <- standardize_data(carsData[,2:11])
Ycars <- carsData[,1]

# generate matrices for each predictor for storing values of all models
cylEst <- matrix(0,200,3, byrow = TRUE)
colnames(cylEst) <- c("full ls", "lasso", "ridge")
dispEst <- matrix(0,200,3, byrow = TRUE)
colnames(dispEst) <- c("full ls", "lasso", "ridge")
hpEst <- matrix(0,200,3, byrow = TRUE)
colnames(hpEst) <- c("full ls", "lasso", "ridge")
dratEst <- matrix(0,200,3, byrow = TRUE)
colnames(dratEst) <- c("full ls", "lasso", "ridge")
wtEst <- matrix(0,200,3, byrow = TRUE)
colnames(wtEst) <- c("full ls", "lasso", "ridge")
qsecEst <- matrix(0,200,3, byrow = TRUE)
colnames(qsecEst) <- c("full ls", "lasso", "ridge")
vsEst <- matrix(0,200,3, byrow = TRUE)
colnames(vsEst) <- c("full ls", "lasso", "ridge")
amEst <- matrix(0,200,3, byrow = TRUE)
colnames(amEst) <- c("full ls", "lasso", "ridge")
gearEst <- matrix(0,200,3, byrow = TRUE)
colnames(gearEst) <- c("full ls", "lasso", "ridge")
carbEst <- matrix(0,200,3, byrow = TRUE)
colnames(carbEst) <- c("full ls", "lasso", "ridge")

for(iter in 1:200){
  idx <- sample(ncar, (2*ncar)/3,replace = TRUE)
  trainData <- Xcars[idx,]
  trainResp <- Ycars[idx]
  
  olsData <- as.data.frame(cbind(trainData,trainResp))
  names(olsData)[ncol(olsData)] <- c("mpg")
  olsEst <- summary(lm(mpg~.-1, olsData))$coefficients[,'Estimate']
  t <- 0.44*sum(abs(olsEst))
  lassoEst <- lassoImp(trainData, as.matrix(trainResp), t)
  ridgeEst <-lm.ridge(mpg~.-1,data=olsData, lambda=0.44)$coef
  
  cylEst[iter,]   <- matrix(c(olsEst[1],lassoEst[1],ridgeEst[1]),3,1)
  dispEst[iter,]  <- matrix(c(olsEst[2],lassoEst[2],ridgeEst[2]),3,1)
  hpEst[iter,]    <- matrix(c(olsEst[3],lassoEst[3],ridgeEst[3]),3,1)
  dratEst[iter,]  <- matrix(c(olsEst[4],lassoEst[4],ridgeEst[4]),3,1)
  wtEst[iter,]    <- matrix(c(olsEst[5],lassoEst[5],ridgeEst[5]),3,1)
  qsecEst[iter,]  <- matrix(c(olsEst[6],lassoEst[6],ridgeEst[6]),3,1)
  vsEst[iter,]    <- matrix(c(olsEst[7],lassoEst[7],ridgeEst[7]),3,1)
  amEst[iter,]   <- matrix(c(olsEst[8],lassoEst[8],ridgeEst[8]),3,1)
  gearEst[iter,]  <- matrix(c(olsEst[9],lassoEst[9],ridgeEst[9]),3,1)
  carbEst[iter,]  <- matrix(c(olsEst[10],lassoEst[10],ridgeEst[10]),3,1)
}


# plots.
par(mfcol=c(3,1))

boxplot(cylEst, main ="box plot of Cylinders")
boxplot(dispEst, main ="box plot of Displacement")
boxplot(hpEst, main ="box plot of HorsePower")
boxplot(dratEst, main ="box plot of Rear Axle Ratio")
boxplot(wtEst, main ="box plot of Weight")
boxplot(qsecEst, main ="box plot of qsec")
boxplot(vsEst, main ="box plot of V/S")
boxplot(amEst, main ="box plot of am")
boxplot(gearEst, main ="box plot of Gears")
boxplot(carbEst, main ="box plot of Carburetors")




