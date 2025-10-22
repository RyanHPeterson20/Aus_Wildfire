#figures from models and model validation

#notes: 
## create figures associated with models

#TODO: create sections for figures (e.g. numbered)
#1. 
# . Adj. R^2 
# . K-fold CV plots
# . RMSE plots (see line ..)
# . CRPS plots
#n. Interval Score Plots


#libraries
suppressMessages( library(grid)) #gridlines between plots
suppressMessages( library(scales)) #for adjusting opacity

suppressMessages( library(RColorBrewer)) #colorRampPalette
suppressMessages( library( fields)) #for set.panel() and others 

#data/model imports
setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("validation_refits.rda") #refits and validation (BIC)
load("validation_refitsEBIC.rda") #refits and validation (eBIC)
load("validation_kfold.rda") #kfold cv for both BIC and eBIC

## --- setup --- ##

#TODO: move model fits up here, if needed (check later)

#TODO: get the changes in the number of terms for eBIC models



## --- Main --- ##

# . (n.) Predictions and Intervals




# . (n.) adj R^2 (in-sample)
#TODO: finalize the below figures

#Get all seasons
#NE Aus 
NE.cons.adjR2 <- NULL
NE.vary.adjR2 <- NULL
for (j in 1:3) {
  cons.adjR2 <- numeric(19)
  vary.adjR2 <- numeric(19)
  for (i in 1:19) {
    cons.lm.temp <- NErefit.new[[2]][[i]]
    vary.lm.temp <- NErefit.new[[3]][[i]]

    cons.adjR2[i] <- summary(cons.lm.temp[[j]])$adj.r.squared
    vary.adjR2[i] <- summary(vary.lm.temp[[j]])$adj.r.squared
  }
  NE.cons.adjR2[[j]] <- cons.adjR2
  NE.vary.adjR2[[j]] <- vary.adjR2
}

#NE Aus group 1
plot(1:19, NE.cons.adjR2[[1]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[1]], pch = 17, col = "forestgreen")

#NE Aus group 2
plot(1:19, NE.cons.adjR2[[2]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[2]], pch = 17, col = "forestgreen")

#NE Aus group 3
plot(1:19, NE.cons.adjR2[[3]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[3]], pch = 17, col = "forestgreen")

#SE Aus 
SE.cons.adjR2 <- NULL
SE.vary.adjR2 <- NULL
for (j in 1:3) {
  cons.adjR2 <- numeric(19)
  vary.adjR2 <- numeric(19)
  for (i in 1:19) {
    cons.lm.temp <- SErefit.new[[2]][[i]]
    vary.lm.temp <- SErefit.new[[3]][[i]]
    
    cons.adjR2[i] <- summary(cons.lm.temp[[j]])$adj.r.squared
    vary.adjR2[i] <- summary(vary.lm.temp[[j]])$adj.r.squared
  }
  SE.cons.adjR2[[j]] <- cons.adjR2
  SE.vary.adjR2[[j]] <- vary.adjR2
}

#SE Aus group 1
plot(1:19, SE.cons.adjR2[[1]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[1]], pch = 17, col = "forestgreen")

#SE Aus group 2
plot(1:19, SE.cons.adjR2[[2]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[2]], pch = 17, col = "forestgreen")

#SE Aus group 3
plot(1:19, SE.cons.adjR2[[3]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[3]], pch = 17, col = "forestgreen")


## 2019-2020 Season only
#NE Aus
#Group 1
NE.cons.adjR2[[1]][19]
NE.vary.adjR2[[1]][19]
#Group 2
NE.cons.adjR2[[2]][19]
NE.vary.adjR2[[2]][19]
#Group 3
NE.cons.adjR2[[3]][19]
NE.vary.adjR2[[3]][19]

#SE Aus
#Group 1
SE.cons.adjR2[[1]][19]
SE.vary.adjR2[[1]][19]
#Group 2
SE.cons.adjR2[[2]][19]
SE.vary.adjR2[[2]][19]
#Group 3
SE.cons.adjR2[[3]][19]
SE.vary.adjR2[[3]][19]


# . (n.) K-Fold CV
#TODO: add legend to each plot, then create output (.png)

#NE Aus
NE1.kcv <- lapply(NE.kcv.ebic, function(x) x[[1]])
NE2.kcv <- lapply(NE.kcv.ebic, function(x) x[[2]])
NE3.kcv <- lapply(NE.kcv.ebic, function(x) x[[3]])

NE.kcv.range <- range(NE1.kcv, NE2.kcv, NE3.kcv)
gamma.lab <- sapply(1:12, function(x) paste0("Gamma: ", round(gamma.seq[x], 3)))

plot(1:12, NE1.kcv, type = "b", pch = 16, col = "firebrick3",
     ylab = "K-fold CV", xlab = "",
     ylim = NE.kcv.range, axes = FALSE)
lines(1:12, NE2.kcv, type = "b", pch =17, col = "royalblue4")
lines(1:12, NE3.kcv, type = "b", pch =17, col = "forestgreen")
box()
axis(1, labels = gamma.lab, at = 1:12, cex.axis = 0.67)
axis(2)  
title("NE Aus K-fold CV")

#SE Aus
SE1.kcv <- lapply(SE.kcv.ebic, function(x) x[[1]])
SE2.kcv <- lapply(SE.kcv.ebic, function(x) x[[2]])
SE3.kcv <- lapply(SE.kcv.ebic, function(x) x[[3]])

SE.kcv.range <- range(SE1.kcv, SE2.kcv, SE3.kcv)
gamma.lab <- sapply(1:12, function(x) paste0("Gamma: ", round(gamma.seq[x], 3)))

plot(1:12, SE1.kcv, type = "b", pch = 16, col = "firebrick3",
     ylab = "K-fold CV", xlab = "",
     ylim = SE.kcv.range, axes = FALSE)
lines(1:12, SE2.kcv, type = "b", pch =17, col = "royalblue4")
lines(1:12, SE3.kcv, type = "b", pch =17, col = "forestgreen")
box()
axis(1, labels = gamma.lab, at = 1:12, cex.axis = 0.67)
axis(2)  
title("SE Aus K-fold CV")




# . (number.) RMSE plots
## RMSE setup
NErmse <- NEvalid[[1]]
NErmse.eBIC <- NEvalid.eBIC[[1]]

## -- base : model fits

#NE Aus group 1
NEgroup1.base.rmse <- unlist(lapply(NErmse, function(x) x[[1]][1]))

#NEbase.rmse.eBIC <- matrix(NA, ncol = 19)
k <- 3
NEbase.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[1]][1]))
k <- 9
NEbase.rmse.eBIC2 <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[1]][1]))

#boxplot
rmse.range <- range(NEgroup1.base.rmse, NEbase.rmse.eBIC, NEbase.rmse.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))


boxplot(NEgroup1.base.rmse, NEbase.rmse.eBIC, NEbase.rmse.eBIC2,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Base Model : RMSE", adj = 0)


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.rmse.eBIC2 <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
rmse.range <- range(NEgroup1.const.rmse, NEgroup1.const.rmse.eBIC, NEgroup1.const.rmse.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.rmse, NEgroup1.const.rmse.eBIC, NEgroup1.const.rmse.eBIC2,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : RMSE", adj = 0)


#NE Aus Group 2
NEgroup2.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
rmse.range <- range(NEgroup2.const.rmse, NEgroup2.const.rmse.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.rmse, NEgroup2.const.rmse.eBIC,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : RMSE", adj = 0)


#NE Aus Group 3
NEgroup3.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][3]))

#eBIC setup
k <- 8
NEgroup3.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
rmse.range <- range(NEgroup3.const.rmse, NEgroup3.const.rmse.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.rmse, NEgroup3.const.rmse.eBIC,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : RMSE", adj = 0)



# . (number.) CRPS plots
## CPRS setup
NEcprs <- NEvalid[[2]]
NEcprs.eBIC <- NEvalid.eBIC[[2]]


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.cprs.eBIC2 <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
cprs.range <- range(NEgroup1.const.cprs, NEgroup1.const.cprs.eBIC, NEgroup1.const.cprs.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.cprs, NEgroup1.const.cprs.eBIC, NEgroup1.const.cprs.eBIC2,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : CRPS", adj = 0)


#NE Aus Group 2
NEgroup2.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
cprs.range <- range(NEgroup2.const.cprs, NEgroup2.const.cprs.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.cprs, NEgroup2.const.cprs.eBIC,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : CRPS", adj = 0)


#NE Aus Group 3
NEgroup3.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][3]))

#eBIC setup
k <- 8
NEgroup3.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
cprs.range <- range(NEgroup3.const.cprs, NEgroup3.const.cprs.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.cprs, NEgroup3.const.cprs.eBIC,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : CRPS", adj = 0)





# . (number.) Interval Score plots
## ints setup
NEints <- NEvalid[[3]]
NEints.eBIC <- NEvalid.eBIC[[3]]


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.ints <- unlist(lapply(NEints, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.ints.eBIC2 <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
ints.range <- range(NEgroup1.const.ints, NEgroup1.const.ints.eBIC, NEgroup1.const.ints.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.ints, NEgroup1.const.ints.eBIC, NEgroup1.const.ints.eBIC2,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : Interval Score", adj = 0)


#NE Aus Group 2
NEgroup2.const.ints <- unlist(lapply(NEints, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
ints.range <- range(NEgroup2.const.ints, NEgroup2.const.ints.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.ints, NEgroup2.const.ints.eBIC,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : Interval Score", adj = 0)


#NE Aus Group 3
NEgroup3.const.ints <- unlist(lapply(NEints, function(x) x[[2]][3]))
#eBIC setup
k <- 8
NEgroup3.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
ints.range <- range(NEgroup3.const.ints, NEgroup3.const.ints.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.ints, NEgroup3.const.ints.eBIC,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : Interval Score", adj = 0)


