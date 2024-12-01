# for the placement of functions used in lasso_cv.R

#prediction and PI

predict.fusedlasso <- function(object, lambda, y_new, X_new) {
  
  
  y <- object$y
  X <- object$X
  beta <- coef(test_object, lambda = lambda)$beta
  
  y_full <- X %*% beta
  
  sigma_hat <- sum((y-y_full)^2)/(n-1)
  
  n <- length(y)
  p <- length(which(beta != 0))
  
  y_new <- as.numeric(y_new) #change to as.numeric()
  X_new <- as.matrix(X_new)
  
  y_hat <- X_new %*% beta
  
  
  #note (X^t X)^{-1} is singular
  A <- t(X) %*% X
  
  X_pseudo <- ginv(A)
  
  
  #quick test on 
  #x_1 <- X_new[1, ]
  #x_2 <- X_new[2, ]
  #x_3 <- X_new[3, ]
  
  h_full <- X_new %*% X_pseudo %*% t(X_new)
  #h_1 <- t(x_1) %*% X_pseudo %*% x_1
  #h_2 <- t(x_2) %*% X_pseudo %*% x_2
  #h_3 <- t(x_3) %*% X_pseudo %*% x_3
  
  #assuming normality
  group_margin <- 1.96*sqrt(sigma_hat*(1 + diag(h_full)))
  temp_upper <- y_hat + group_margin
  temp_lower <- y_hat - group_margin
  
  pred_df <- data.frame(y_new, y.hat = y_hat, PI.upper = temp_upper, 
                        PI.lower = temp_lower)
  colnames(pred_df) <- c("y"," y.hat", "PI.upper", "PI.lower")
  rownames(pred_df) <- NULL
  
  return(pred_df)
}




#R2





#