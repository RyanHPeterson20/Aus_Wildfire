
#following the .ncl method for cmip5 model data

#TODO: wip for pIOD index functions for SST anoms and PCA/EOF


#following the .ncl method for cmip5 model data
#currently only provides SON anoms
sst.anoms <- function(A, son.only = FALSE, quadratic = FALSE){
  
  A.new <- aperm(A, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]
  
  #get dimensions from reordered array
  nx <- dim(A.new)[3]
  ny <- dim(A.new)[2] 
  nt <- dim(A.new)[1]
  
  #get SON mean (get 1 sst per year (as son mean))
  if (son.only) {
    #if only provided with 3 months (son) per year
    nyears <- nt %/% 3
    A.temp <- array(A.new, dim = c(3, nyears, ny, nx))
    A.son <- apply(A.temp[1:3, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset
  }else {
    nyears <- nt %/% 12
    #reshape data array
    A.temp <- array(A.new, dim = c(12, nyears, ny, nx))
    A.son <- apply(A.temp[9:11, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset
  }
  
  rm(A.temp)
  
  #update time dim
  nt <- dim(A.son)[1]
  
  #remove mean
  var <- matrix(A.son, nrow = nt, ncol = ny * nx) #from earlier work on detrend
  
  #get full period mean
  var.mean <- colMeans(var) #TODO: visualize this spatially (transform to [lon, lat])
  var <- sweep(var, 2, colMeans(var), "-") #TODO: check how this might hold up with irregular spatial masks
  
  #de-trend
  #TODO: figure out how to work with irregular NAs in time/space (some years have different spatial NAs)
  test.na <- which(!is.na(var), arr.ind = TRUE) 
  true.index <- unique(test.na[,2])
  
  A.resid <- matrix(NA, nrow = nt, ncol = ny*nx)
  A.coef <- matrix(NA, nrow = 2, ncol = ny*nx)
  
  var.x <- seq_len(nt) #time "covariates"
  var.y <- var[ ,true.index] #response variable (sst anoms)
  
  #TODO: add in an elif for linear and quadratic, with and without intercept
  
  if (quadratic) {
    var.fit <- lm(var.y ~ var.x + I(var.x^2))
    A.coef[,true.index] <- var.fit$coefficients[1:2, ]
  } else {
    var.fit <- lm(var.y ~ var.x)
    A.coef[,true.index] <- var.fit$coefficients
  }
  
  A.resid[,true.index] <- var.fit$residuals
  #TODO: update this handle quadratic fit
  #A.coef[,true.index] <- var.fit$coefficients
  
  return(list(mean = var.mean,
              anom = A.resid,
              coef = A.coef))
}


sst.anoms.new <- function(A, season = "son", quadratic = FALSE){
  
  A.new <- aperm(A, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]
  
  #get dimensions from reordered array
  nx <- dim(A.new)[3]
  ny <- dim(A.new)[2] 
  nt <- dim(A.new)[1]
  
  #get seasonal mean (get 1 sst per year (as seasonal mean))  
  nyears <- nt %/% 12
  #reshape data array
  A.temp <- array(A.new, dim = c(12, nyears, ny, nx))
  
  #TODO: setup season input check, this is pretty ad-hoc and needs some rails
  
  if (season == "djf") {
    #djf setup
    A.season <- array(NA, dim = c(nyears - 1, ny, nx))
    
    for (y in 1:(nyears - 1)) {
      dec <- A.temp[12, y, , ]       # December of year y
      jan <- A.temp[1, y + 1, , ]    # January of year y+1
      feb <- A.temp[2, y + 1, , ]    # February of year y+1
      # Compute mean across the three months
      A.season[y, , ] <- apply(abind(dec, jan, feb, along = 0), c(2, 3), mean, na.rm = TRUE)
    }
  } else {
    season.mon <- switch (season,
                          son = 9:11, jja = 6:8, mam = 3:5)
    
    A.season <- apply(A.temp[season.mon, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset  
    
  }
  
  rm(A.temp)
  
  #update time dim
  nt <- dim(A.season)[1]
  
  #remove mean
  var <- matrix(A.season, nrow = nt, ncol = ny * nx) #from earlier work on detrend
  
  #get full period mean
  var.mean <- colMeans(var) #TODO: visualize this spatially (transform to [lon, lat])
  var <- sweep(var, 2, colMeans(var), "-") #TODO: check how this might hold up with irregular spatial masks
  
  #de-trend
  #TODO: figure out how to work with irregular NAs in time/space (some years have different spatial NAs)
  test.na <- which(!is.na(var), arr.ind = TRUE) 
  true.index <- unique(test.na[,2])
  
  A.resid <- matrix(NA, nrow = nt, ncol = ny*nx)
  A.coef <- matrix(NA, nrow = 2, ncol = ny*nx)
  
  var.x <- seq_len(nt) #time "covariates"
  var.y <- var[ ,true.index] #response variable (sst anoms)
  
  #TODO: add in an elif for linear and quadratic, with and without intercept
  
  if (quadratic) {
    var.fit <- lm(var.y ~ var.x + I(var.x^2))
    A.coef[,true.index] <- var.fit$coefficients[1:2, ]
  } else {
    var.fit <- lm(var.y ~ var.x)
    A.coef[,true.index] <- var.fit$coefficients
  }
  
  A.resid[,true.index] <- var.fit$residuals
  #TODO: update this handle quadratic fit
  #A.coef[,true.index] <- var.fit$coefficients
  
  return(list(mean = var.mean,
              anom = A.resid,
              coef = A.coef))
}



#pca/eof analysis function:
##sst.anom as [time, lat, lon]
sst.eof <- function(sst.anom, kmode){

  nt <- dim(sst.anom)[1]
  ny <- dim(sst.anom)[2]
  nx <- dim(sst.anom)[3]
  
  X <- matrix(sst.anom, nrow = nt, ncol = ny * nx)
  
  #get rid of NA's (masked locs)
  keep <- colSums(is.finite(X)) == nt
  X.new  <- X[, keep]
  
  svd.temp <- svd(X.new)
  
  #svd outputs
  U <- svd.temp$u[ ,1:kmode]
  D <- svd.temp$d[1:kmode]
  V <- svd.temp$v[ ,1:kmode]
  
  #outputs
  EOF.temp <- V #EOF spatial pattern
  PC.temp <- U %*% diag(D)
  per.temp <- D^2 / sum(svd.temp$d^2) 
  
  #finalize the spatial eof (pca)
  V_eof <- matrix(NA, nrow = ny * nx, ncol = kmode)
  
  #add in lsm sea data
  V_eof[keep, ] <- EOF.temp
  
  return(list(EOF = V_eof,
              PC = PC.temp,
              percent = per.temp))
}


#son averages
## sst as [lon, lat, time]
son.avg <- function(sst, son.only = FALSE){
  ##sst as [lon, lat, time]
  A.new <- aperm(sst, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]
  
  #get dimensions from reordered array
  nx <- dim(A.new)[3]
  ny <- dim(A.new)[2] 
  nt <- dim(A.new)[1]
  
  #get SON mean (get 1 sst per year (as son mean))
  if (son.only) {
    #if only provided with 3 months (son) per year
    nyears <- nt %/% 3
    A.temp <- array(A.new, dim = c(3, nyears, ny, nx))
    A.son <- apply(A.temp[1:3, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset
  }else {
    nyears <- nt %/% 12
    #reshape data array
    A.temp <- array(A.new, dim = c(12, nyears, ny, nx))
    A.son <- apply(A.temp[9:11, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset
  }
  
  return(aperm(A.son, perm = c(3, 2, 1))) #return as [lon, lat, time]
}

