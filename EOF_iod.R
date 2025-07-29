
#EOF (PCA) function

#TODO: finalize params and outputs

#parameters
## A: weekly field data (not anomaly) in [lon, lat, time] (or [nx, ny, nt])
## kmode: number of modes from PCA
## amode: selects for the anomaly mode, TRUE = for single year climatology, FALSE (else) = all climate (full dataset)
## lsm: land sea mask, 1 = sea (valid), 0 = land (invalid) as [nx, ny] ([lon, lat])
## weekly: true = for weekly data, false = for monthly data

#output:
# V_eof is spatial pattern (as num matrix)
# pc_ts is temporal pattern (as num matrix)
# per is proportion for each mode (as num vector)


eof_iod <- function(A, kmode, amode = TRUE, lsm, weekly = TRUE){
  
  Q_net <- A #as [nx, ny, nt] ([lon, lat, time]) 
  
  Q_net_1 <- aperm(Q_net, perm = c(2, 1, 3))  # reorder as [ny, nx, nt] ([lat, lon, time])
  
  dims <- dim(Q_net_1) #get dimensions from reordered array
  ny <- dims[1] 
  nx <- dims[2]
  nt <- dims[3]
  
  # Overwrite Q_net and clear Q_net_1
  Q_net <- Q_net_1
  rm(Q_net_1)
  
  # 'flatten' to [ny*nx, nt]
  np <- nx * ny
  var <- matrix(Q_net, nrow = np, ncol = nt)
  
  #anomaly (yearly)
  if (amode) {
    if (weekly) {
      week_mean <- matrix(NA, nrow = np, ncol = 52)
      for (j in 1:52) {
        # get all indices for week j (j in 1-52)
        cols_j <- seq(j, nt, by = 52)
        
        # get weekly mean for each spatial point 
        week_mean[ ,j] <- rowMeans(var[, cols_j, drop = FALSE])
        
        # number of time points for j
        ncols <- length(cols_j)
        
        # subtract weekly mean
        var[, cols_j] <- var[, cols_j] - matrix(week_mean[, j], nrow = np, ncol = ncols)
      }
      
    } else{ #monthly
      month_mean <- matrix(NA, nrow = np, ncol = 12)
      
      for (k in 1:12) {
        # get all indices for week k (k in 1-12)
        cols_k <- seq(k, nt, by = 12)
        
        # get monthly mean for each spatial point 
        month_mean[ ,k] <- rowMeans(var[, cols_k, drop = FALSE])
        
        # number of time points for k
        ncols <- length(cols_j)  
        
        # subtract monthly mean
        var[, cols_k] <- var[, cols_k] - matrix(week_mean[, k], nrow = np, ncol = ncols)
      }
    }
  } else{ #full dataset anom
    var <- sweep(var, 1, rowMeans(var), "-")
  }
  
  #land-sea-mask (lsm)
  lsm_vec <- as.vector(lsm) #flatten mask from [nx, ny] to [np]  
  
  # keep sea points (where mask == 1)
  is_sea <- lsm_vec == 1
  var_eof <- var[is_sea, ] #select for only sea surface temp

  
  #find PCA, using prcomp (centered data)
  pca <- prcomp(t(var_eof), center = TRUE, scale. = TRUE )

  #EOF spatial [np, kmod]
  pc_EOF <- pca$rotation[, 1:kmode]
  
  #pca time series
  pc_ts <- pca$x[, 1:kmode]
  
  #percent/proportion for each mode
  per <- pca$sdev[1:kmode]^2 / sum(pca$sdev^2)
  
  #finalize the spatial eof (pca)
  V_eof <- matrix(NA, nrow = np, ncol = kmode)
  
  #add in lsm sea data
  V_eof[lsm_vec == 1, ] <- pc_EOF
  
  # TODO: update for better names
  return(list(
    V_eof = V_eof,
    pc_ts = matrix2,
    per = per
  )) 
}
