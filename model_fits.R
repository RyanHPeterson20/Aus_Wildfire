
#notes:
##-next steps in finalizing the model 
##-add duplicate to finalproject_code folder

#TODO: add more libraries as needed
#library 
suppressMessages( library(glmnet)) #test ridge regression for coefs

suppressMessages( library( Metrics)) #measurement metrics

#data import
#TODO: collect data into fewer locations
setwd("~/CO_AUS/Aus_CO-main/Interactions")
load( "Data/ne_data.rda")
load( "Data/se_data.rda")
load( "Data/bounded_data.rda")
load( "Data/data_matrix.rda")
load( "Data/lag_list.rda")
load( "Data/data_quantile.rda")

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("IOD_lag.rda") #IOD data

#function load

