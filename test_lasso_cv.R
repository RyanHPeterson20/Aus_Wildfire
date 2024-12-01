

#current method is to see how we perform CV with gen/fused lasso. 
#reproducing cv.trendfilter but for regression
setwd("~/CO_AUS/Aus_CO-main")
source("genlasso_func.R")

#Update with trendfilter

#folds 
k <- 5

#data and 
y <- object_test$y
X <- object_test$X
D = Matrix(D_new, sparse=TRUE)

#direct parameters passed from fusedLasso
maxsteps=2000
minlam=0
rtol=1e-7
btol=1e-7
eps=1e-4

n <- length(y)

foldid <- c(0,rep(Seq(1,k),n-2)[Seq(1,n-2)], 0) #fold id

lambda <- object_test$lambda

#not needed (keep for now)
ord = object_test$ord #NULL
ord = 0
pos = object_test$pos #NULL

#0 matrix for each fold
cvall <- matrix(0,k,length(lambda))


i <- 1
#for (i in genlasso::Seq(1,k)){
otr <- which(foldid != i) #training id 
ntr <- length(otr) #length of training data
ytr <- y[otr]
Xtr <- X[otr, ]


#example: function call
##Dtr = getDtfPosSparse(ntr,ord,ptr)
##out = dualpathWideSparse(ytr,Dtr,NULL,approx,Inf,min(lambda),rtol,btol,verbose)
## (Have to do this manually now) since The parent funtion will return the proper beta, fit, y, bls
#- from above we would find beta, fit, y, bls via trendfilter/fusedlasso
##out$beta = as.matrix(ytr - t(Dtr)%*%out$u)

#dualpathFusedX(y,X,D,approx,maxsteps,minlam,rtol,btol,eps,verbose)
object_test$gamma

out_new <- fusedlasso(y = ytr, X = Xtr, D_new, gamma = 0.85, minlam = min(lambda))


length(out_new$lambda) #note: smaller than the original object

#trying to get this to work (we need to find the correct: t(out_new$pathobjs$D2))
#below is not needed
#as.matrix(ytr - t(out_new$pathobjs$D2)%*%out_new$u)
#out_fusedLasso <- genlasso::dualpathFusedX(ytr, Xtr, D, approx = FALSE, maxsteps, minlam, rtol, btol, eps)


lambda_new <- out_new$lambda

##fits
b_new <- coef(out_new, lambda = lambda)$beta

ote <- which(foldid == i) #testing id
yte <- matrix(y[ote], length(ote), length(lambda)) #testing data
Xte <- X[ote, ]
#Xte <- matrix(X[ote, ], length(ote), length(lambda))

#probably not needed
ilo <- which((Seq(1,n)%in%(ote-1))[otr])
ihi <- which((Seq(1,n)%in%(ote+1))[otr])

test_pred <- Xte%*%b_new
colMeans((yte-test_pred)^2)

##now let's construct a better cv functions

#TODO: once this is working move over to lasso_valid_functions.R
#parameters:
##-object (fused lasso object for full fit)
##-k (number of folds)
##-

cv.fusedlasso <- function(object, k =5, D){
  cl = match.call()
  
  y <- object$y
  X <- object$X
  D_new <- Matrix(D, sparse=TRUE)
  n <- length(y)
  
  foldid <- c(0, rep(Seq(1, k), n-2)[Seq(1, n-2)], 0) #fold id
  lambda <- object$lambda
  
  cvall <- matrix(0, nrow = k, ncol = length(lambda))
  
  for (i in Seq(1, k)) {
    
    otr <- which(foldid != i) #training id 
    ntr <- length(otr) #length of training data
    ytr <- y[otr] #training response
    Xtr <- X[otr, ] #training predictors
    
    gamma <- object$gamma
    
    out_new <- fusedlasso(y = ytr, X = Xtr, D, gamma = gamma, minlam = min(lambda))
    
    b_new <- coef(out_new, lambda = lambda)$beta
    
    ote <- which(foldid == i) #testing id
    yte <- matrix(y[ote], length(ote), length(lambda)) #testing data
    Xte <- X[ote, ]
    
    pred_new <- Xte%*%b_new
    cvall[i,] <- colMeans((yte-pred_new)^2) #cv mse
    
  }
  
  #directly from cv.trendfilter.R
  cverr <- colMeans(cvall)
  cvse <- apply(cvall,2,sd)/sqrt(k)
  
  names(cverr) = names(cvse) = round(lambda,3)
  i0 = which.min(cverr)
  lam.min = lambda[i0]
  lam.1se = max(lambda[cverr<=cverr[i0]+cvse[i0]])
  i.min = which(lambda==lam.min)
  i.1se = which(lambda==lam.1se)
  
  out = list(err=cverr,se=cvse,mode="lambda",lambda=lambda,
             lambda.min=lam.min,lambda.1se=lam.1se,i.min=i.min,i.1se=i.1se,call=cl)
  
  class(out) = c("cv.trendfilter", "list")  
  return(out)
}


#test above function
out_cv <- cv.fusedlasso(object_test, k = 5, D = D_new)
out_cv$lambda.min
plot(out_cv)



#test new cv/validation functions

#current parameters include:
##object (As test_object)
##lambda (As test_lambda)
##y_new (As test_resp)
##X_new (As test_preds)

y <- test_object$y
X <- test_object$X
beta <- coef(test_object, lambda = test_lambda)$beta

y_full <- X %*% beta

MSE <- mean((y-y_full)^2)

n <- length(y)
p <- length(which(beta != 0))

y_new <- test_resp #change to as.numeric()
X_new <- as.matrix(test_preds)

y_hat <- X_new %*% beta


#note (X^t X)^{-1} is singular
A <- t(X) %*% X

X_pseudo <- ginv(A)


#quick test on 
x_1 <- X_new[1, ]
x_2 <- X_new[2, ]
x_3 <- X_new[3, ]

h_full <- X_new %*% X_pseudo %*% t(X_new)
h_1 <- t(x_1) %*% X_pseudo %*% x_1
h_2 <- t(x_2) %*% X_pseudo %*% x_2
h_3 <- t(x_3) %*% X_pseudo %*% x_3

SSR_hat <- sum((y_new - y_hat)^2)
