#helper functions for model refits

#refit ramp function
#parameters
## test.ramp: ramp model object
## X_1: predictor data object as matrix

refit_ramp <- function(test.ramp, X_1){
  
  main.terms <- test.ramp$mainInd
  main.terms <- colnames(X_1)[main.terms]
  
  interactions <- test.ramp$interInd
  if (!is.null(interactions)){
    for (i in 1:length(interactions)){
      
      this.term <- interactions[i]
      these.subterms <- as.integer(strsplit(this.term, "X")[[1]][2:3])
      these.subterms <- colnames(X_1)[these.subterms]
      
      if (these.subterms[1] != these.subterms[2]){
        this.term <- paste(these.subterms, collapse = ":")
      } else {
        this.term <- paste0("I(", these.subterms[1], "^2)")
      }
      
      interactions[i] <- this.term
    }
  }
  
  #refit with lm
  if (!is.null(interactions)){
    model.string <- paste( paste(main.terms, collapse = " + "),
                           paste(interactions, collapse = " + "),
                           sep = " + ")
  } else {
    model.string <- paste(main.terms, collapse = " + ")
  }
  
  model.string <- paste0("co ~ ", model.string)
  
  return(model.string)
}



terms_only <- function(test.ramp, X_1){
  main.terms <- test.ramp$mainInd
  main.terms <- colnames(X_1)[main.terms]
  
  interactions <- test.ramp$interInd
  if (!is.null(interactions)){
    for (i in 1:length(interactions)){
      
      this.term <- interactions[i]
      these.subterms <- as.integer(strsplit(this.term, "X")[[1]][2:3])
      these.subterms <- colnames(X_1)[these.subterms]
      
      if (these.subterms[1] != these.subterms[2]){
        this.term <- paste(these.subterms, collapse = ":")
      } else {
        this.term <- paste0("I(", these.subterms[1], "^2)")
      }
      
      interactions[i] <- this.term
    }
  }
  return(c(main.terms, interactions))
}
