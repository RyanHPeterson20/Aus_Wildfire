
#R script for finding SST anomalies

#libraries
#.nc files
suppressMessages(library(ncdf4))
suppressMessages(library(terra))

# date mgmt
suppressMessages(library(lubridate))
suppressMessages(library(abind))

#visualization and interpolation
suppressMessages(library(fields))
suppressMessages(library(colorspace))

#matlab function (e.g. detrend)
suppressMessages(library(pracma))


#load in data
setwd("~/CO_AUS/Aus_CO-main/pIOD")
load("sst_base.rda")

#setup 
lsm.IOD.array <- array(lsm.IOD, dim = dim(sst.OISST.new))

sst.OISST.masked <- ifelse(lsm.IOD.array == 1, sst.OISST.new, NA)

#TODO: check for NA's across the entire period (that is, find a uniform LSM)
#check on SODA
##which(is.na(sst.interp.SODA))
#Check on GODAS

#internal functions


##----Main----##

#following the .ncl method for cmip5 model data
#TODO: functionalize if this works
A <- sst.OISST.masked #as [lon, lat, time]

A.new <- aperm(A, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]

dims <- dim(A.new) #get dimensions from reordered array
nx <- dims[3]
ny <- dims[2] 
nt <- dims[1]

#get SON mean
nyears <- nt %/% 12
#reshape data array
A.temp <- array(A.new, dim = c(12, nyears, ny, nx))
A.son <- apply(A.temp[9:11, , , ], c(2,3,4), mean, na.rm = TRUE) #already at correct subset
rm(A.temp)

#update time dim
nt <- dim(A.son)[1]

#remove mean
var <- matrix(A.son, nrow = nt, ncol = ny * nx) #from earlier work on detrend
#get full period mean
var.mean <- colMeans(var) #TODO: visualize this spatially (transform to [lon, lat])

var <- sweep(var, 2, colMeans(var), "-") #TODO: check how this might hold up with irregular spatial masks

#detrend
#TODO: test linear detrend and plot slope
test.na <- which(!is.na(var), arr.ind = TRUE)

A.resid <- matrix(NA, nrow = nt, ncol = ny*nx)
A.coef <- matrix(NA, nrow = 2, ncol = ny*nx)

#TODO: figure out how to work with irregular NAs in time/space (some years have different spatial NAs)
true.index <- unique(test.na[,2])

var.x <- seq_len(nt) #time "covariates"
var.y <- var[ ,true.index] #response variable (sst anoms)

var.fit <- lm(var.y ~ var.x)

A.resid[,true.index] <- var.fit$residuals
A.coef[,true.index] <- var.fit$coefficients





#Visualizations
## strong years 1994, 1997, 2006 (index: 13, 16, 25)
## moderate years 1982, 1987, 2015 (index: 1, )
## along with spatial mean, linear coeffs

#array reformatting code
#TODO: delete when done, (or move to function) 
temp <- array(var, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.diff <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]


#Spatial Mean
#using var.mean
spat.mean <- array(var.mean, dim = c(ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]

image.plot(list(x = lon.values, y = rev(lat.values), z = t(spat.mean)), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)


#Spatial detrend coefs
temp <- array(A.coef, dim = c(2, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.coefs <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,2]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)



#1993 SON average
son.avg <- aperm(A.son, c(3, 2, 1))

image.plot(list(x = lon.values, y = rev(lat.values), z = son.avg[,,13]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")#, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)

#1993 SON anom

temp <- array(var, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.diff <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.diff[,,13]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")#, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)

#1993 SON anom detrended

temp <- array(A.resid, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.resid <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,13]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")#, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)





##------TEST SECTION-----##
#test viz (remove or move down when done):
#transform from [time, ny*nx] to [time, lat, lon]
temp <- array(var, dim = c(nt, ny, nx))  # t(var) is [nt, ny*nx]
spat.diff <- aperm(temp, c(3, 2, 1))  # [ny, nx, nt]
##spat.out <- aperm(Q_net_1, c(2, 1, 3)) #  Reverse the earlier aperm to get back to [nx, ny, nt]
#transform from [time, lat, lon] to [lon, lat, time]
spat.out <- aperm(A.son, c(3, 2 ,1))


#set.panel(1,3)
image.plot(list(x = lon.values, y = rev(lat.values), z = spat.diff[,,16]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")#, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)


image.plot(list(x = lon.values, y = rev(lat.values), z = sst.OISST.masked[,,23]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat", main = paste0("OISST: ", format(as.Date(times[23+1]), "%B %Y")))
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = sst.interp.SODA[,,25]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat", main = paste0("SODA3.3.1: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)


dev.off()
