# functions for grouping of response and predictor variables:


#grouped response set-up functions
##TODO: refine these functions to handle any response clusters.

#NE resp vectors
NEresp_grouping <- function(NEAus_mat, j = (1:19)){
  
  NEAus_mat <- lapply(NEAus_mat, function(mat) mat[j, ])
  
  NEAus1_resp <- as.vector(NEAus_mat[[1]])
  NEAus2_resp <- as.vector(NEAus_mat[[2]])
  NEAus3_resp <- as.vector(NEAus_mat[[3]])
  NEAus4_resp <- as.vector(NEAus_mat[[4]])
  NEAus5_resp <- as.vector(NEAus_mat[[5]])
  NEAus6_resp <- as.vector(NEAus_mat[[6]])
  
  NEAus_vec <- list(NEAus1_resp, NEAus2_resp, NEAus3_resp,
                    NEAus4_resp, NEAus5_resp, NEAus6_resp)
  
  return(NEAus_vec)
}

#SE resp vectors
SEresp_grouping <- function(SEAus_mat, j = (1:19)){
  
  SEAus_mat <- lapply(SEAus_mat, function(mat) mat[j, ])
  
  SEAus1_resp <- as.vector(SEAus_mat[[1]])
  SEAus2_resp <- as.vector(SEAus_mat[[2]])
  SEAus3_resp <- as.vector(SEAus_mat[[3]])
  SEAus4_resp <- as.vector(SEAus_mat[[4]])
  SEAus5_resp <- as.vector(SEAus_mat[[5]])
  
  SEAus_vec <- list(SEAus1_resp, SEAus2_resp, SEAus3_resp,
                    SEAus4_resp, SEAus5_resp)
  
  return(SEAus_vec) 
}  


#lag predictor set-up
#TODO: fine tune these functions for different groupings/clustering
#goal is to enter in the grouping and get the appropriate outputs
#-currently manually set-up groups

## NE lag predictors
NElag_grouping <- function(NE_laglist, j = (1:19)){
  #row inclusion
  #j <- (1:19)
  
  
  #NE group 1: 35-37
  NE1_lag <- NE_laglist$`Week  35`[j, -(1:2)]
  NE1_lag <- rbind(NE1_lag, NE_laglist$`Week  36`[j, -(1:2)], NE_laglist$`Week  37`[j, -(1:2)])
  
  #NE group 2: 38-42
  NE2_lag <- NE_laglist$`Week  38`[j,-(1:2)]
  NE2_lag <- rbind(NE2_lag, NE_laglist$`Week  39`[j,-(1:2)], NE_laglist$`Week  40`[j,-(1:2)],
                   NE_laglist$`Week  41`[j,-(1:2)], NE_laglist$`Week  42`[j,-(1:2)])
  
  #NE group 3: 43-46
  NE3_lag <- NE_laglist$`Week  43`[j,-(1:2)]
  NE3_lag <- rbind(NE3_lag, NE_laglist$`Week  44`[j,-(1:2)], NE_laglist$`Week  45`[j,-(1:2)],
                   NE_laglist$`Week  46`[j,-(1:2)])
  
  #NE group 4: 47-51
  NE4_lag <- NE_laglist$`Week  47`[j,-(1:2)]
  NE4_lag <- rbind(NE4_lag, NE_laglist$`Week  48`[j,-(1:2)], NE_laglist$`Week  49`[j,-(1:2)],
                   NE_laglist$`Week  50`[j,-(1:2)], NE_laglist$`Week  51`[j,-(1:2)])
  
  #NE group 5: 52, 1-3
  NE5_lag <- NE_laglist$`Week  52`[j,-(1:2)]
  NE5_lag <- rbind(NE5_lag, NE_laglist$`Week  1`[j,-(1:2)], NE_laglist$`Week  2`[j,-(1:2)],
                   NE_laglist$`Week  3`[j,-(1:2)])
  
  #NE group 6: 4-14
  NE6_lag <- NE_laglist$`Week  4`[j,-(1:2)]
  NE6_lag <- rbind(NE6_lag, NE_laglist$`Week  5`[j,-(1:2)], NE_laglist$`Week  6`[j,-(1:2)],
                   NE_laglist$`Week  7`[j,-(1:2)], NE_laglist$`Week  8`[j,-(1:2)],
                   NE_laglist$`Week  9`[j,-(1:2)], NE_laglist$`Week  10`[j,-(1:2)],
                   NE_laglist$`Week  11`[j,-(1:2)], NE_laglist$`Week  12`[j,-(1:2)],
                   NE_laglist$`Week  13`[j,-(1:2)], NE_laglist$`Week  14`[j,-(1:2)])
  
  NEAus_preds <- list(NE1_lag, NE2_lag, NE3_lag, NE4_lag, NE5_lag, NE6_lag)
  return(NEAus_preds)
}  

## SE lag predictors
SElag_grouping <- function(SE_laglist, j = (1:19)){
  
  #SE group 1: 35-37
  SE1_lag <- SE_laglist$`Week  35`[j,-(1:2)]
  SE1_lag <- rbind(SE1_lag, SE_laglist$`Week  36`[j,-(1:2)], SE_laglist$`Week  37`[j,-(1:2)])
  
  #SE group 2: 38-41
  SE2_lag <- SE_laglist$`Week  38`[j,-(1:2)]
  SE2_lag <- rbind(SE2_lag, SE_laglist$`Week  39`[j,-(1:2)], SE_laglist$`Week  40`[j,-(1:2)],
                   SE_laglist$`Week  41`[j,-(1:2)])
  
  #SE group 3: 42-50
  SE3_lag <- SE_laglist$`Week  42`[j,-(1:2)]
  SE3_lag <- rbind(SE3_lag, SE_laglist$`Week  43`[j,-(1:2)], SE_laglist$`Week  44`[j,-(1:2)],
                   SE_laglist$`Week  45`[j,-(1:2)], SE_laglist$`Week  46`[j,-(1:2)],
                   SE_laglist$`Week  47`[j,-(1:2)], SE_laglist$`Week  48`[j,-(1:2)],
                   SE_laglist$`Week  49`[j,-(1:2)], SE_laglist$`Week  50`[j,-(1:2)])
  
  #SE group 4: 51-52,1-2
  SE4_lag <- SE_laglist$`Week  51`[j,-(1:2)]
  SE4_lag <- rbind(SE4_lag, SE_laglist$`Week  52`[j,-(1:2)], SE_laglist$`Week  1`[j,-(1:2)],
                   SE_laglist$`Week  2`[j,-(1:2)])
  
  #SE group 5: 3-14
  SE5_lag <- SE_laglist$`Week  3`[j,-(1:2)]
  SE5_lag <- rbind(SE5_lag, SE_laglist$`Week  4`[j,-(1:2)], SE_laglist$`Week  5`[j,-(1:2)],
                   SE_laglist$`Week  6`[j,-(1:2)], SE_laglist$`Week  7`[j,-(1:2)],
                   SE_laglist$`Week  8`[j,-(1:2)], SE_laglist$`Week  9`[j,-(1:2)],
                   SE_laglist$`Week  10`[j,-(1:2)], SE_laglist$`Week  11`[j,-(1:2)],
                   SE_laglist$`Week  12`[j,-(1:2)], SE_laglist$`Week  13`[j,-(1:2)],
                   SE_laglist$`Week  14`[j,-(1:2)])
  
  
  
  SEAus_preds <- list(SE1_lag, SE2_lag, SE3_lag, SE4_lag, SE5_lag)
  return(SEAus_preds)
}  