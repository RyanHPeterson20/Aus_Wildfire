# for the placement of functions used in lasso_cv.R

#prediction and PI

predict.fusedlasso <- function(object, test_year, lambda) {
  
  y <- test_object$y
  X <- test_object$X
  beta <- coef(test_object, lambda = test_lambda)$beta
  
  n <- length(y)
  p <- length(which(beta != 0))
}




#R2





#