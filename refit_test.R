##

#testing function for getting refits/BIC

#base coefs (from hierNet.path)
get_coefs <- function(path_group, max_index, resp, preds, preds_quant){
  
  #extract coefs from path object
  #TODO: set object (list) for coefs
  term_list <- list()
  for (k in 1:max_index) {
    #k <- 11
    #TODO: redo this part as its own function so I can call 
    temp_bp <- path_group$bp[,k]
    temp_bn <- path_group$bn[,k]
    
    temp_th <- path_group$th[,,k]
    
    #internals
    y = as.numeric(resp)
    X = cbind(as.matrix(preds[ ,1:260]),
              as.matrix(preds_quant[ ,1:104])  )
    
    
    coef_names <- c(colnames(preds), paste0( "IND_", colnames(preds_quant[ ,1:104]) ) )
    colnames(X) <- coef_names
    
    #main effects
    main_effect <- temp_bp - temp_bn #gets main effect coefs
    main_terms <- which(main_effect != 0, arr.ind = TRUE)
    mains <- colnames(X)[main_terms]
    
    main_df <- data.frame(mains, round(main_effect[main_terms], 5))
    
    rownames(main_df) <- NULL
    colnames(main_df) <- c("Main_Effect", "Coef")
    
    #interaction effects
    interact_effect <- temp_th
    interact_terms <- which(interact_effect != 0, arr.ind = TRUE)
    
    #TODO: add elif for interactions
    if (nrow(interact_terms) != 0) {
      interact_names <- matrix(coef_names[interact_terms], nrow = nrow(interact_terms), ncol = ncol(interact_terms))
      
      interactions <-  paste0(interact_names[,2], ":", interact_names[,1])
      for (i in 1:nrow(interact_names)) {
        interact_subterms <- interact_names[i, ]
        if (interact_subterms[1] != interact_subterms[2]) {
          #condition for interactions
          this_term <- paste(interact_subterms, collapse = ":")
        } else {
          #condition for squared terms
          this_term <- paste0("I(", interact_subterms[1], "^2)")
        }
        
        interactions[i] <- this_term
      }
      interact_df <- data.frame(interactions, round(interact_effect[interact_terms], 5) )
      
      rownames(interact_df) <- NULL
      colnames(interact_df) <- c("Interact_Effect", "Coef")
      #update main + interaction effects
      term_list[[k]] <- list(main_df, interact_df)
      
    } else {
      #only main effects update
      term_list[[k]] <- list(main_df)
    }
  }
  return(term_list)
}



#assign max lambda index at the start
#refit with BIC function
refit_bic <- function(path_group, max_index, ebic.gamma, resp, preds, preds_quant){


  #lambda_mse <- which(cv_group$lamlist == cv_group$lamhat)
  #lambda_1se <- which(cv_group$lamlist == cv_group$lamhat.1se)
  
  refit_lm <- list() #refit lm
  refit_lmcoef <- list() #lm coefs
  refit_ridge <- list() #refit ridge
  refit_ridgecoef <- list() #ridge coefs
  refit_ridgecv <- list() #ridge cv (for R^2)
  BIC_matrix <- matrix(NA, ncol = 5)
  colnames(BIC_matrix) <- c("lambda", "lm.BIC", "ridge.BIC", "lm.eBIC", "ridge.eBIC")
  
  
  for (k in 1:max_index) {
    #k <- 11
    #TODO: redo this part as its own function so I can call 
    temp_bp <- path_group$bp[,k]
    temp_bn <- path_group$bn[,k]
    
    temp_th <- path_group$th[,,k]
    
  
    
    #internals
    y = as.numeric(resp)
    X = cbind(as.matrix(preds[ ,1:260]),
              as.matrix(preds_quant[ ,1:104])  )
    
    
    coef_names <- c(colnames(preds), paste0( "IND_", colnames(preds_quant[ ,1:104]) ) )
    colnames(X) <- coef_names
    
    #main effects
    main_effect <- temp_bp - temp_bn
    main_terms <- which(main_effect != 0, arr.ind = TRUE)
    mains <- colnames(X)[main_terms]
    
    #mains
    #main_effect[main_terms]
    
    #length(main_terms)
    #length(interactions)
    
    #interaction effects
    interact_effect <- temp_th
    interact_terms <- which(interact_effect != 0, arr.ind = TRUE)
    
    #interact_effect[interact_terms]
    
    #TODO: add elif for interactions
    if (nrow(interact_terms) != 0) {
      interact_names <- matrix(coef_names[interact_terms], nrow = nrow(interact_terms), ncol = ncol(interact_terms))
      
      interactions <-  paste0(interact_names[,2], ":", interact_names[,1])
      for (i in 1:nrow(interact_names)) {
        interact_subterms <- interact_names[i, ]
        if (interact_subterms[1] != interact_subterms[2]) {
          #condition for interactions
          this_term <- paste(interact_subterms, collapse = ":")
        } else {
          #condition for squared terms
          this_term <- paste0("I(", interact_subterms[1], "^2)")
        }
        
        interactions[i] <- this_term
      }
      
    }
  
    #setup model for other fits
    if (nrow(interact_terms) != 0) {
      model_string <- paste( paste(mains, collapse = " + "),
                             paste(interactions, collapse = " + "),
                             sep = " + ")
    } else  {
      model_string <- paste(mains, collapse = " + ")
    }
    
    
    model_string <- paste0("y ~ ", model_string)
    
    
    #elif
    #linear model
    data_df <- as.data.frame(cbind(y,X))
    lm_fit <- lm(formula(model_string), data = data_df)
    lm_coef <- coef(lm_fit)
    
    #ridge model
    X_df <- as.data.frame(X)
    f <- as.formula(model_string)
    X_new <- model.matrix(f, X_df)
    
    #TODO: update below with elif for lambda.min and lambda.1se (include in )
    set.seed(300)
    ridge_cv <- cv.glmnet(X_new, y, alpha = 0.00, nfolds = 5)
    ridge_fit <- glmnet(X_new, y, alpha = 0.00)
    coef_pred <- predict(ridge_fit, s = ridge_cv$lambda.1se, type = "coefficients")
    par_vec <- coef_pred@i + 1 
    ridge_coef <- coef_pred[par_vec, ] #coefs from ridge
    
    lm_resid <- y - predict(lm_fit, X_df)
    ridge_resid <- y - predict(ridge_fit, X_new, s = ridge_cv$lambda.1se, type = "response")
    
    #BIC
    lm_BIC <- length(y)*log(mean(lm_resid^2)) + log(length(y))*length(lm_coef[-1])
    ridge_BIC <- length(y)*log(mean(ridge_resid^2)) + log(length(y))*length(ridge_coef[-1])
    
    #eBIC 
    #get p effective.
    p <- dim(X)[2]
    df.main <- length(main_terms)
    p_eff <- p + df.main *(df.main + 1)/2
    
    lm_eBIC <- lm_BIC + 2 * ebic.gamma * log(choose(p_eff, length(lm_coef[-1])))
    ridge_eBIC <- ridge_BIC +  2 * ebic.gamma * log(choose(p_eff, length(ridge_coef[-1])))
    
    refit_lm[[paste0("LamIndex_", k)]] <- lm_fit
    refit_lmcoef[[paste0("LamIndex_", k)]] <- lm_coef
    refit_ridge[[paste0("LamIndex_", k)]] <- ridge_fit
    refit_ridgecoef[[paste0("LamIndex_", k)]] <- ridge_coef
    refit_ridgecv[[paste0("LamIndex_", k)]] <- ridge_cv
    
    BIC_matrix <- rbind(BIC_matrix, c(path_group$lamlist[k], lm_BIC, ridge_BIC, lm_eBIC, ridge_eBIC))
  }
  
  BIC_matrix <- BIC_matrix[-1, ]
  
  return_list <- list(refit_lm, refit_lmcoef, refit_ridge, refit_ridgecoef, BIC_matrix, refit_ridgecv)
  
  return(return_list)
}


