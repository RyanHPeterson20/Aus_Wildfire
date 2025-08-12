
#R script for finding SST anomalies

#libraries
#.nc files
library(ncdf4)
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
nrow(test.na)
true.index <- unique(test.na[,2])
var.y <- var[ ,true.index]



#transform from [time, ny*nx] to [time, lat, lon]
temp <- array(var, dim = c(nt, ny, nx))  # t(var) is [nt, ny*nx]
spat.diff <- aperm(temp, c(3, 2, 1))  # [ny, nx, nt]
##spat.out <- aperm(Q_net_1, c(2, 1, 3)) #  Reverse the earlier aperm to get back to [nx, ny, nt]
#transform from [time, lat, lon] to [lon, lat, time]
spat.out <- aperm(A.son, c(3, 2 ,1))

#test viz (remove or move down when done):


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
