##internal functions for the genlasso 
#note: this is only for internal use and not for publication


Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}


#TODO: expand for determining R^2 and other metrics.
#OR copy some of this code to produce more interesting statistics.

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