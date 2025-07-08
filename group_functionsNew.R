

NEresp_grouping <- function(NEAus_mat, j = (1:19)){
  
  NEAus_mat <- lapply(NEAus_mat, function(mat) mat[j, ])
  
  NEAus1_resp <- as.vector(NEAus_mat[[1]])
  NEAus2_resp <- as.vector(NEAus_mat[[2]])

  NEAus_vec <- list(NEAus1_resp, NEAus2_resp)
  
  return(NEAus_vec)
}


#SE resp vectors
SEresp_grouping <- function(SEAus_mat, j = (1:19)){
  
  SEAus_mat <- lapply(SEAus_mat, function(mat) mat[j, ])
  
  SEAus1_resp <- as.vector(SEAus_mat[[1]])
  SEAus2_resp <- as.vector(SEAus_mat[[2]])
  #SEAus3_resp <- as.vector(SEAus_mat[[3]])
  
  SEAus_vec <- list(SEAus1_resp, SEAus2_resp)
  
  return(SEAus_vec) 
}  



## NE lag predictors
NElag_grouping <- function(NE_laglist, j = (1:19)){
  #row inclusion
  #j <- (1:19)
  
  #NE group 1: Weeks 35-46 Index (1-12) 
  NE1_lag <- NE_laglist$`Week  35`[j, -(1:2)]
  NE1_lag <- rbind(NE1_lag, NE_laglist$`Week  36`[j, -(1:2)], NE_laglist$`Week  37`[j, -(1:2)],
                   NE_laglist$`Week  38`[j,-(1:2)], NE_laglist$`Week  39`[j,-(1:2)], NE_laglist$`Week  40`[j,-(1:2)],
                   NE_laglist$`Week  41`[j,-(1:2)], NE_laglist$`Week  42`[j,-(1:2)],
                   NE_laglist$`Week  43`[j,-(1:2)],
                   NE_laglist$`Week  44`[j,-(1:2)], NE_laglist$`Week  45`[j,-(1:2)],
                   NE_laglist$`Week  46`[j,-(1:2)])
                   
  #NE group 2  Weeks 47-52; 1-14 Index (13-32)
  NE2_lag <- NE_laglist$`Week  47`[j,-(1:2)]
  NE2_lag <- rbind( NE2_lag, NE_laglist$`Week  48`[j,-(1:2)], NE_laglist$`Week  49`[j,-(1:2)],
                    NE_laglist$`Week  50`[j,-(1:2)], NE_laglist$`Week  51`[j,-(1:2)], 
                    NE_laglist$`Week  52`[j,-(1:2)], NE_laglist$`Week  1`[j,-(1:2)], 
                    NE_laglist$`Week  2`[j,-(1:2)], NE_laglist$`Week  3`[j,-(1:2)], 
                    NE_laglist$`Week  4`[j,-(1:2)], NE_laglist$`Week  5`[j,-(1:2)], 
                    NE_laglist$`Week  6`[j,-(1:2)], NE_laglist$`Week  7`[j,-(1:2)], 
                    NE_laglist$`Week  8`[j,-(1:2)], NE_laglist$`Week  9`[j,-(1:2)], 
                    NE_laglist$`Week  10`[j,-(1:2)], NE_laglist$`Week  11`[j,-(1:2)], 
                    NE_laglist$`Week  12`[j,-(1:2)], NE_laglist$`Week  13`[j,-(1:2)], 
                    NE_laglist$`Week  14`[j,-(1:2)])
  

  NEAus_preds <- list(NE1_lag, NE2_lag)
  return(NEAus_preds)
} 


## SE lag predictors
SElag_grouping <- function(SE_laglist, j = (1:19)){

  #SE group 1:  Weeks 35-42 Index (1-8)
  SE1_lag <- SE_laglist$`Week  35`[j,-(1:2)]
  SE1_lag <- rbind(SE1_lag, SE_laglist$`Week  36`[j,-(1:2)], SE_laglist$`Week  37`[j,-(1:2)], 
                   SE_laglist$`Week  38`[j,-(1:2)], SE_laglist$`Week  39`[j,-(1:2)], 
                   SE_laglist$`Week  40`[j,-(1:2)], SE_laglist$`Week  41`[j,-(1:2)],
                   SE_laglist$`Week  42`[j,-(1:2)])
                   
                   
                   

  #SE group 2: Weeks 43-52; 1-14 Index (9-32)
  SE2_lag <- SE_laglist$`Week  43`[j,-(1:2)]
  SE2_lag <- rbind(SE2_lag, SE_laglist$`Week  44`[j,-(1:2)],
                   SE_laglist$`Week  45`[j,-(1:2)], SE_laglist$`Week  46`[j,-(1:2)],
                   SE_laglist$`Week  47`[j,-(1:2)], SE_laglist$`Week  48`[j,-(1:2)],
                   SE_laglist$`Week  49`[j,-(1:2)], SE_laglist$`Week  50`[j,-(1:2)], 
                   SE_laglist$`Week  51`[j,-(1:2)], SE_laglist$`Week  52`[j,-(1:2)], 
                   SE_laglist$`Week  1`[j,-(1:2)], SE_laglist$`Week  2`[j,-(1:2)],
                   SE_laglist$`Week  3`[j,-(1:2)], SE_laglist$`Week  4`[j,-(1:2)], SE_laglist$`Week  5`[j,-(1:2)],
                   SE_laglist$`Week  6`[j,-(1:2)], SE_laglist$`Week  7`[j,-(1:2)],
                   SE_laglist$`Week  8`[j,-(1:2)], SE_laglist$`Week  9`[j,-(1:2)],
                   SE_laglist$`Week  10`[j,-(1:2)], SE_laglist$`Week  11`[j,-(1:2)],
                   SE_laglist$`Week  12`[j,-(1:2)], SE_laglist$`Week  13`[j,-(1:2)],
                   SE_laglist$`Week  14`[j,-(1:2)])
  
  

  
  SEAus_preds <- list(SE1_lag, SE2_lag)
  return(SEAus_preds)
} 



#Lag groupings for ward.D2 (sq euclidean)

NEresp_3group <- function(NEAus_mat, j = (1:19)){
  
  NEAus_mat <- lapply(NEAus_mat, function(mat) mat[j, ])
  
  NEAus1_resp <- as.vector(NEAus_mat[[1]])
  NEAus2_resp <- as.vector(NEAus_mat[[2]])
  NEAus3_resp <- as.vector(NEAus_mat[[3]])
  
  NEAus_vec <- list(NEAus1_resp, NEAus2_resp, NEAus3_resp)
  
  return(NEAus_vec)
}


#SE resp vectors
SEresp_3group <- function(SEAus_mat, j = (1:19)){
  
  SEAus_mat <- lapply(SEAus_mat, function(mat) mat[j, ])
  
  SEAus1_resp <- as.vector(SEAus_mat[[1]])
  SEAus2_resp <- as.vector(SEAus_mat[[2]])
  SEAus3_resp <- as.vector(SEAus_mat[[3]])
  
  SEAus_vec <- list(SEAus1_resp, SEAus2_resp, SEAus3_resp)
  
  return(SEAus_vec) 
}  


#NE lag predictors
NElag_old <- function(NE_laglist, j = (1:19)){
  #row inclusion
  #j <- (1:19)
  
  #NE group 1: Weeks 35-46 Index (1-12)
  NE1_lag <- NE_laglist$`Week  35`[j, -(1:2)]
  NE1_lag <- rbind(NE1_lag, NE_laglist$`Week  36`[j, -(1:2)], NE_laglist$`Week  37`[j, -(1:2)],
                   NE_laglist$`Week  38`[j,-(1:2)], NE_laglist$`Week  39`[j,-(1:2)], NE_laglist$`Week  40`[j,-(1:2)],
                   NE_laglist$`Week  41`[j,-(1:2)], NE_laglist$`Week  42`[j,-(1:2)],
                   NE_laglist$`Week  43`[j,-(1:2)],
                   NE_laglist$`Week  44`[j,-(1:2)], NE_laglist$`Week  45`[j,-(1:2)],
                   NE_laglist$`Week  46`[j,-(1:2)])
  
  #NE group 2  Weeks 47-51 Index (13-17)
  NE2_lag <- NE_laglist$`Week  47`[j,-(1:2)]
  NE2_lag <- rbind( NE2_lag, NE_laglist$`Week  48`[j,-(1:2)], NE_laglist$`Week  49`[j,-(1:2)],
                    NE_laglist$`Week  50`[j,-(1:2)], NE_laglist$`Week  51`[j,-(1:2)]) 
  
  #NE group 3 Week 52; 1-14 Index (18-32)
  NE3_lag <- NE_laglist$`Week  52`[j,-(1:2)]
  NE3_lag <- rbind(NE3_lag, NE_laglist$`Week  1`[j,-(1:2)], 
                   NE_laglist$`Week  2`[j,-(1:2)], NE_laglist$`Week  3`[j,-(1:2)], 
                   NE_laglist$`Week  4`[j,-(1:2)], NE_laglist$`Week  5`[j,-(1:2)], 
                   NE_laglist$`Week  6`[j,-(1:2)], NE_laglist$`Week  7`[j,-(1:2)], 
                   NE_laglist$`Week  8`[j,-(1:2)], NE_laglist$`Week  9`[j,-(1:2)], 
                   NE_laglist$`Week  10`[j,-(1:2)], NE_laglist$`Week  11`[j,-(1:2)],
                   NE_laglist$`Week  12`[j,-(1:2)], NE_laglist$`Week  13`[j,-(1:2)], 
                   NE_laglist$`Week  14`[j,-(1:2)])
  
  NEAus_preds <- list(NE1_lag, NE2_lag, NE3_lag)
  return(NEAus_preds)
}

## SE lag predictors
SElag_old <- function(SE_laglist, j = (1:19)){
  
  #SE group 1:  Weeks 35-50 Index (1-16)
  SE1_lag <- SE_laglist$`Week  35`[j,-(1:2)]
  SE1_lag <- rbind(SE1_lag, SE_laglist$`Week  36`[j,-(1:2)], SE_laglist$`Week  37`[j,-(1:2)],
                   SE_laglist$`Week  38`[j,-(1:2)], SE_laglist$`Week  39`[j,-(1:2)],
                   SE_laglist$`Week  40`[j,-(1:2)], SE_laglist$`Week  41`[j,-(1:2)],
                   SE_laglist$`Week  42`[j,-(1:2)], SE_laglist$`Week  43`[j,-(1:2)],
                   SE_laglist$`Week  44`[j,-(1:2)],
                   SE_laglist$`Week  45`[j,-(1:2)], SE_laglist$`Week  46`[j,-(1:2)],
                   SE_laglist$`Week  47`[j,-(1:2)], SE_laglist$`Week  48`[j,-(1:2)],
                   SE_laglist$`Week  49`[j,-(1:2)], SE_laglist$`Week  50`[j,-(1:2)])
  

  #SE group 2: Weeks 51-52; 1-2 Index (17-20)
  SE2_lag <-  SE_laglist$`Week  51`[j,-(1:2)]
  SE2_lag <- rbind(SE2_lag, SE_laglist$`Week  52`[j,-(1:2)], 
                   SE_laglist$`Week  1`[j,-(1:2)], SE_laglist$`Week  2`[j,-(1:2)])
  
  
  #SE group 3: Weeks 3-14 Index (21-32)
  SE3_lag <-  SE_laglist$`Week  3`[j,-(1:2)]
  SE3_lag <- rbind(SE3_lag, SE_laglist$`Week  4`[j,-(1:2)], SE_laglist$`Week  5`[j,-(1:2)],
                   SE_laglist$`Week  6`[j,-(1:2)], SE_laglist$`Week  7`[j,-(1:2)],
                   SE_laglist$`Week  8`[j,-(1:2)], SE_laglist$`Week  9`[j,-(1:2)],
                   SE_laglist$`Week  10`[j,-(1:2)], SE_laglist$`Week  11`[j,-(1:2)],
                   SE_laglist$`Week  12`[j,-(1:2)], SE_laglist$`Week  13`[j,-(1:2)],
                   SE_laglist$`Week  14`[j,-(1:2)])
  
  
  SEAus_preds <- list(SE1_lag, SE2_lag, SE3_lag)
  return(SEAus_preds)
} 


## NE lag predictors
NElag_3group <- function(NE_laglist, j = (1:19)){
  #row inclusion
  #j <- (1:19)
  
  #Update group 1 (Weeks 38-46)
  NE1_lag <- NE_laglist$`Week  38`[j, -(1:2)]
  NE1_lag <- rbind(NE1_lag, NE_laglist$`Week  39`[j,-(1:2)], NE_laglist$`Week  40`[j,-(1:2)],
                   NE_laglist$`Week  41`[j,-(1:2)], NE_laglist$`Week  42`[j,-(1:2)],
                   NE_laglist$`Week  43`[j,-(1:2)],
                   NE_laglist$`Week  44`[j,-(1:2)], NE_laglist$`Week  45`[j,-(1:2)],
                   NE_laglist$`Week  46`[j,-(1:2)])
  
  #NE group 2  Weeks 47-51 Index (13-17)
  NE2_lag <- NE_laglist$`Week  47`[j,-(1:2)]
  NE2_lag <- rbind( NE2_lag, NE_laglist$`Week  48`[j,-(1:2)], NE_laglist$`Week  49`[j,-(1:2)],
                    NE_laglist$`Week  50`[j,-(1:2)], NE_laglist$`Week  51`[j,-(1:2)]) 
                    
  #NE group 3 Week 52; 1-14 Index (18-32)
  NE3_lag <- NE_laglist$`Week  52`[j,-(1:2)]
  NE3_lag <- rbind(NE3_lag, NE_laglist$`Week  1`[j,-(1:2)], 
                   NE_laglist$`Week  2`[j,-(1:2)], NE_laglist$`Week  3`[j,-(1:2)], 
                   NE_laglist$`Week  4`[j,-(1:2)], NE_laglist$`Week  5`[j,-(1:2)], 
                   NE_laglist$`Week  6`[j,-(1:2)], NE_laglist$`Week  7`[j,-(1:2)], 
                   NE_laglist$`Week  8`[j,-(1:2)], NE_laglist$`Week  9`[j,-(1:2)], 
                   NE_laglist$`Week  10`[j,-(1:2)], NE_laglist$`Week  11`[j,-(1:2)],
                   NE_laglist$`Week  12`[j,-(1:2)], NE_laglist$`Week  13`[j,-(1:2)], 
                   NE_laglist$`Week  14`[j,-(1:2)])
  
  NEAus_preds <- list(NE1_lag, NE2_lag, NE3_lag)
  return(NEAus_preds)
} 





## SE lag predictors
SElag_3group <- function(SE_laglist, j = (1:19)){
  
  #SE group 1:  Weeks 38-50 Index (1-16)
  SE1_lag <- SE_laglist$`Week  38`[j,-(1:2)]
  SE1_lag <- rbind(SE1_lag, SE_laglist$`Week  39`[j,-(1:2)],
                   SE_laglist$`Week  40`[j,-(1:2)], SE_laglist$`Week  41`[j,-(1:2)],
                   SE_laglist$`Week  42`[j,-(1:2)], SE_laglist$`Week  43`[j,-(1:2)],
                   SE_laglist$`Week  44`[j,-(1:2)],
                   SE_laglist$`Week  45`[j,-(1:2)], SE_laglist$`Week  46`[j,-(1:2)],
                   SE_laglist$`Week  47`[j,-(1:2)], SE_laglist$`Week  48`[j,-(1:2)],
                   SE_laglist$`Week  49`[j,-(1:2)], SE_laglist$`Week  50`[j,-(1:2)])

  
  #SE group 2: Weeks 51-52; 1-2 Index (17-20)
  SE2_lag <-  SE_laglist$`Week  51`[j,-(1:2)]
  SE2_lag <- rbind(SE2_lag, SE_laglist$`Week  52`[j,-(1:2)], 
                   SE_laglist$`Week  1`[j,-(1:2)], SE_laglist$`Week  2`[j,-(1:2)])
                   

  #SE group 3: Weeks 3-14 Index (21-32)
  SE3_lag <-  SE_laglist$`Week  3`[j,-(1:2)]
  SE3_lag <- rbind(SE3_lag, SE_laglist$`Week  4`[j,-(1:2)], SE_laglist$`Week  5`[j,-(1:2)],
                   SE_laglist$`Week  6`[j,-(1:2)], SE_laglist$`Week  7`[j,-(1:2)],
                   SE_laglist$`Week  8`[j,-(1:2)], SE_laglist$`Week  9`[j,-(1:2)],
                   SE_laglist$`Week  10`[j,-(1:2)], SE_laglist$`Week  11`[j,-(1:2)],
                   SE_laglist$`Week  12`[j,-(1:2)], SE_laglist$`Week  13`[j,-(1:2)],
                   SE_laglist$`Week  14`[j,-(1:2)])
  
  
  SEAus_preds <- list(SE1_lag, SE2_lag, SE3_lag)
  return(SEAus_preds)
} 



#without boundaries (both removal)
## NE lag predictors
NElag_new <- function(NE_laglist, j = (1:19)){
  #row inclusion
  #j <- (1:19)
  
  #Update group 1 (Weeks 38-46)
  NE1_lag <- NE_laglist$`Week  38`[j, -(1:2)]
  NE1_lag <- rbind(NE1_lag, NE_laglist$`Week  39`[j,-(1:2)], NE_laglist$`Week  40`[j,-(1:2)],
                   NE_laglist$`Week  41`[j,-(1:2)], NE_laglist$`Week  42`[j,-(1:2)],
                   NE_laglist$`Week  43`[j,-(1:2)],
                   NE_laglist$`Week  44`[j,-(1:2)], NE_laglist$`Week  45`[j,-(1:2)],
                   NE_laglist$`Week  46`[j,-(1:2)])
  
  #NE group 2  Weeks 47-51 Index (13-17)
  NE2_lag <- NE_laglist$`Week  47`[j,-(1:2)]
  NE2_lag <- rbind( NE2_lag, NE_laglist$`Week  48`[j,-(1:2)], NE_laglist$`Week  49`[j,-(1:2)],
                    NE_laglist$`Week  50`[j,-(1:2)], NE_laglist$`Week  51`[j,-(1:2)]) 
  
  #NE group 3 Week 52; 1-14 Index (18-32)
  NE3_lag <- NE_laglist$`Week  52`[j,-(1:2)]
  NE3_lag <- rbind(NE3_lag, NE_laglist$`Week  1`[j,-(1:2)], 
                   NE_laglist$`Week  2`[j,-(1:2)], NE_laglist$`Week  3`[j,-(1:2)], 
                   NE_laglist$`Week  4`[j,-(1:2)], NE_laglist$`Week  5`[j,-(1:2)], 
                   NE_laglist$`Week  6`[j,-(1:2)], NE_laglist$`Week  7`[j,-(1:2)], 
                   NE_laglist$`Week  8`[j,-(1:2)], NE_laglist$`Week  9`[j,-(1:2)], 
                   NE_laglist$`Week  10`[j,-(1:2)], NE_laglist$`Week  11`[j,-(1:2)])
  
  NEAus_preds <- list(NE1_lag, NE2_lag, NE3_lag)
  return(NEAus_preds)
} 


## SE lag predictors
SElag_new <- function(SE_laglist, j = (1:19)){
  
  #SE group 1:  Weeks 38-50 Index (1-16)
  SE1_lag <- SE_laglist$`Week  38`[j,-(1:2)]
  SE1_lag <- rbind(SE1_lag, SE_laglist$`Week  39`[j,-(1:2)],
                   SE_laglist$`Week  40`[j,-(1:2)], SE_laglist$`Week  41`[j,-(1:2)],
                   SE_laglist$`Week  42`[j,-(1:2)], SE_laglist$`Week  43`[j,-(1:2)],
                   SE_laglist$`Week  44`[j,-(1:2)],
                   SE_laglist$`Week  45`[j,-(1:2)], SE_laglist$`Week  46`[j,-(1:2)],
                   SE_laglist$`Week  47`[j,-(1:2)], SE_laglist$`Week  48`[j,-(1:2)],
                   SE_laglist$`Week  49`[j,-(1:2)], SE_laglist$`Week  50`[j,-(1:2)])
  
  
  #SE group 2: Weeks 51-52; 1-2 Index (17-20)
  SE2_lag <-  SE_laglist$`Week  51`[j,-(1:2)]
  SE2_lag <- rbind(SE2_lag, SE_laglist$`Week  52`[j,-(1:2)], 
                   SE_laglist$`Week  1`[j,-(1:2)], SE_laglist$`Week  2`[j,-(1:2)])
  
  
  #SE group 3: Weeks 3-14 Index (21-32)
  SE3_lag <-  SE_laglist$`Week  3`[j,-(1:2)]
  SE3_lag <- rbind(SE3_lag, SE_laglist$`Week  4`[j,-(1:2)], SE_laglist$`Week  5`[j,-(1:2)],
                   SE_laglist$`Week  6`[j,-(1:2)], SE_laglist$`Week  7`[j,-(1:2)],
                   SE_laglist$`Week  8`[j,-(1:2)], SE_laglist$`Week  9`[j,-(1:2)],
                   SE_laglist$`Week  10`[j,-(1:2)], SE_laglist$`Week  11`[j,-(1:2)])
  
  
  SEAus_preds <- list(SE1_lag, SE2_lag, SE3_lag)
  return(SEAus_preds)
} 



#update code to be more general.


