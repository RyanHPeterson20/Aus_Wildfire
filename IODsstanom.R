
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
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales)) #for adjusting opacity

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
#currently only provides SON anoms
sst.anoms <- function(A){
    
  A.new <- aperm(A, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]
  
  #get dimensions from reordered array
  nx <- dim(A.new)[3]
  ny <- dim(A.new)[2] 
  nt <- dim(A.new)[1]
  
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
  
  #de-trend
  #TODO: figure out how to work with irregular NAs in time/space (some years have different spatial NAs)
  test.na <- which(!is.na(var), arr.ind = TRUE) 
  true.index <- unique(test.na[,2])
  
  A.resid <- matrix(NA, nrow = nt, ncol = ny*nx)
  A.coef <- matrix(NA, nrow = 2, ncol = ny*nx)
  
  var.x <- seq_len(nt) #time "covariates"
  var.y <- var[ ,true.index] #response variable (sst anoms)
  
  #TODO: add in an elif for linear and quadratic
  var.fit <- lm(var.y ~ var.x)
  
  A.resid[,true.index] <- var.fit$residuals
  A.coef[,true.index] <- var.fit$coefficients
  
  return(list(mean = var.mean,
              anom = A.resid,
              coef = A.coef))
}

#pca/eof analysis function:
sst.eof <- function(sst.anom, kmode){
  #sst.anom as [time, lat, lon]
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


#sst anomalies SON, for 
##OISST
sst.anom.oisst <- sst.anoms(sst.OISST.masked)
#GODAS
sst.anom.godas <- sst.anoms(sst.array)
#SODA
sst.anom.soda <- sst.anoms(sst.interp.SODA)

#get SON average of SST anomalies
sst.anom.avg <- (sst.anom.oisst$anom + sst.anom.godas$anom + sst.anom.soda$anom)/3

#get EOFs for full Indian Ocean region, then 
#reproduce Fig 1a and 1b: project SST anoms onto EOF1 and EOF2
nt <- 34
nx <- length(lon.values)
ny <- length(lat.values)

sst.anom <- array(sst.anom.avg, dim = c(nt, ny, nx)) 
pca.base <- sst.eof(sst.anom, kmode = 2)

eof.base <- pca.base$EOF
pc.base <- pca.base$PC

pc.base1 <- scale(pc.base[,1], center = TRUE, scale = TRUE)
pc.base2 <- scale(pc.base[,2], center = TRUE, scale = TRUE)

temp <- array(t(eof.base), dim = c(2, ny, nx))  # t(var) is [nt, ny*nx]
#v.eof1 <- temp[1,,]
#v.eof2 <- temp[2,,]
eof.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]

#setup color and range
eof.absmax <- max(abs(eof.spatial), na.rm = TRUE)
eof.range <- c(-eof.absmax, eof.absmax) 

#color setup
cols.ryb <- colorRampPalette(brewer.pal(11, "RdYlBu"))(24)
brks <- seq(eof.range[1], eof.range[2], length.out = length(cols.ryb) + 1)

#TODO: save a final image
image.plot(list(x = lon.values, y = rev(lat.values), z = eof.spatial[,,1]), 
           col = cols.ryb, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = eof.spatial[,,2]), 
           col = cols.ryb, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

#TODO: project SST anoms onto each EOF pattern
eof.mask <- is.finite(eof.base[,1]) #currently looking at eof1

v.eof1 <- eof.base[eof.mask,1] #eof1 vector (may or may not use this)
x.anoms <- sst.anom.avg[,eof.mask] #sst anom data matrix

spatial.proj <- matrix(NA, ncol = 2, nrow = 3600)
proj1 <- t(x.anoms) %*% pc.std.IOD[,1]
proj2 <- t(x.anoms) %*% pc.std.IOD[,2]

spatial.proj[eof.mask, ] <- cbind(proj1, proj2)

temp <- array(t(spatial.proj), dim = c(2, ny, nx))  # t(var) is [nt, ny*nx]
#v.eof1 <- temp[1,,]
#v.eof2 <- temp[2,,]
proj.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]

proj.spatial <- proj.spatial/34

#setup color and range
eof.absmax <- max(abs(proj.spatial), na.rm = TRUE)
eof.range <- c(-eof.absmax, eof.absmax) 

#color setup
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
brks1 <- c(eof.range[1], seq(-0.5, 0.5, 0.1), eof.range[2])
brks <- seq(eof.range[1], eof.range[2], length.out = length(cols.ryb) + 1)

#TODO: save a final image
#TODO: adjust to more "exactly" match
image.plot(list(x = lon.values, y = rev(lat.values), z = proj.spatial[,,1]), 
           col = cols.ryb, breaks = brks1, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = proj.spatial[,,2]), 
           col = cols.ryb, breaks = brks1, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)


#Update for figures 1e and 1f.
#reconstruction (UDV^T) of temp anoms
z1 <- pc.std.IOD[,1] %*% t(proj1)
z2 <- pc.std.IOD[,2] %*% t(proj2)

test.z1 <- colSums(z1)
test.z2 <- colSums(z2)

spatial.anoms <- matrix(NA, ncol = 2, nrow = 3600)
spatial.anoms <- 
spatial.proj[eof.mask, ] <- cbind(proj1, proj2)

temp <- array(t(spatial.proj), dim = c(2, ny, nx))  # t(var) is [nt, ny*nx]
#v.eof1 <- temp[1,,]
#v.eof2 <- temp[2,,]
proj.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]

s.spatial <- (proj.spatial[,,1] + proj.spatial[,,2])
m.spatial <- (proj.spatial[,,1] - proj.spatial[,,2])

#setup color and range
index.absmax <- max(abs(s.spatial), abs(m.spatial), na.rm = TRUE)
index.range <- c(-index.absmax, index.absmax) 

brks.index <- c(index.range[1], seq(-1, 1, 0.2), index.range[2])

image.plot(list(x = lon.values, y = rev(lat.values), z = s.spatial), 
           col = cols.ryb, breaks = brks.index, zlim = index.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = m.spatial), 
           col = cols.ryb, breaks = brks.index, zlim = index.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)



#select for reduced region pIOD 
IOD_maxLon <- 100
IOD_minLon <- 40
IOD_maxLat <- 5
IOD_minLat <- -5

#location range
lon.values.IOD <- lon.values[lon.values <= IOD_maxLon & lon.values >= IOD_minLon]
lon.range.IOD <- range(which(lon.values <= IOD_maxLon & lon.values >= IOD_minLon))

lat.values.IOD <- lat.values[lat.values <= IOD_maxLat & lat.values >= IOD_minLat]
lat.range.IOD <- range(which(lat.values <= IOD_maxLat & lat.values >= IOD_minLat))

sst.anom <- array(sst.anom.avg, dim = c(nt, ny, nx)) 

sst.anom.pIOD <- sst.anom[,lat.range.IOD[1]:lat.range.IOD[2],lon.range.IOD[1]:lon.range.IOD[2]]

#reduced IOD pca
pca.pIOD <- sst.eof(sst.anom.pIOD, kmode = 2)

pc.std.IOD <- scale(pca.pIOD$PC, center = TRUE, scale = TRUE)

pc.std.IOD[,1] <- -pc.std.IOD[,1] 

s.index <- (pc.std.IOD[,1] + pc.std.IOD[,2])/sqrt(2)
m.index <- (pc.std.IOD[,1] - pc.std.IOD[,2])/sqrt(2)

#TODO: save these further down
plot(1:34, s.index, type = "l", col = "firebrick", ylim=c(-2,5))
lines(1:34, m.index, col = "darkgreen")
abline(h=c(1.5, 1.25), lty = 2, col = c("firebrick", "darkgreen"))

plot(pc.std.IOD[,1], pc.std.IOD[,2], pch = 16, xlim = c(-2,4), ylim = c(-3, 5),
     xlab = "PC1", ylab = "PC2", col = "darkblue")
points(pc.std.IOD[c(13,16,25),1], pc.std.IOD[c(13,16,25),2], pch = 16, col = "firebrick")
text(pc.std.IOD[25,1], pc.std.IOD[25,2], "2006", pos = 4, cex = 0.9, col = "firebrick")
text(pc.std.IOD[16,1], pc.std.IOD[16,2], "1997", pos = 3, cex = 0.9, col = "firebrick")
text(pc.std.IOD[13,1], pc.std.IOD[13,2], "1994", pos = 4, cex = 0.9, col = "firebrick")
points(pc.std.IOD[c(1,6,34),1], pc.std.IOD[c(1,6,34),2], pch = 16, col = "darkgreen")
text(pc.std.IOD[34,1], pc.std.IOD[34,2], "2015", pos = 4, cex = 0.9, col = "darkgreen")
text(pc.std.IOD[6,1], pc.std.IOD[6,2], "1987", pos = 4, cex = 0.9, col = "darkgreen")
text(pc.std.IOD[1,1], pc.std.IOD[1,2], "1982", pos = 3, cex = 0.9, col = "darkgreen")
abline(h = 0.5, lty = 2)
abline(v = 1, lty = 2)
abline(v = c(1.1, 1.25), lty = 2, col = "firebrick")


#TODO: move the sst regression to here


#TODO: testing transparent saves
setwd("~/CO_AUS/Aus_CO-main/pIOD/Figures")

#adjusted colors for overlay
png("overlay.png", width = 2000, height = 1200, bg = "transparent", type = "cairo") 
par(mar = c(0,0,0,0), oma = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1:34, s.index, type = "l", col = "orange", ylim=c(-2,5),
     axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE, lwd = 9)
lines(1:34, m.index, col = "lawngreen", lwd = 9)
abline(h=c(1.5, 1.25), lty = 2, col = c("violetred4", "seagreen4"), lwd = 6)
dev.off()

png("PCoverlay_temp.png", width = 1500, height = 1500, bg = "transparent", type = "cairo") 
par(mar = c(0,0,0,0), oma = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(pc.std.IOD[,1], pc.std.IOD[,2], pch = 21, col = "black",
     bg = alpha("cyan",.65), xlim = c(-2,4), ylim = c(-3, 5),
     xlab = "PC1", ylab = "PC2", cex = 5)
points(pc.std.IOD[c(13,16,25),1], pc.std.IOD[c(13,16,25),2], pch = 21, col = "black",
       bg = alpha("red1",.75), cex = 5)
points(pc.std.IOD[c(1,6,34),1], pc.std.IOD[c(1,6,34),2], pch = 21, col = "black",
       bg = alpha("seagreen1",.65), cex = 5)
abline(h = 0.5, lty = 2, lwd = 4)
abline(v = 1, lty = 2, lwd = 4)
abline(v = c(1.1, 1.25), lty = 2, col = "firebrick")
dev.off()

#reproduce Fig 1a and 1b:
#project SST anoms onto EOF1 and EOF2


#
temp <- array(t(V_eof), dim = c(2, ny, nx))  # t(var) is [nt, ny*nx]
#v.eof1 <- -temp[1,,]
#v.eof2 <- temp[2,,]
eof.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]
#eof.spatial <- aperm(Q_net_1, c(2, 1, 3))

#setup color and range
eof.absmax <- max(abs(eof.spatial), na.rm = TRUE)
eof.range <- c(-eof.absmax, eof.absmax) 
roma_col <- rev(divergingx_hcl(n = 48, palette = "RdYlBu"))

brks <- seq(eof.range[1], eof.range[2], length.out = length(roma_col) + 1)

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = -eof.spatial[,,1]), 
           col = roma_col, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)


image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = eof.spatial[,,2]), 
           col = roma_col, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)


#Visualizations
##s-, m-index
years <- 1982:2015
son.dates <- make_date(years, 1, 1) 

setwd("~/CO_AUS/Aus_CO-main/pIOD/Figures")

#pIOD indices
png(filename = "pIODindex_new.png", width = 2600, height = 1300, res = 250)
par(mar = c(5, 5, 3, 5))
plot(son.dates, s.index, type = "l", col = "firebrick3", ylim = c(-2, 5),
     xlab = "Year", ylab = "Index", lwd = 2, xaxt = "n")
axis.Date(1, at = son.dates[seq(4,34, by = 5)], format = "%Y")
#axis.Date(1, at = seq(min(son.dates), max(son.dates), by = "5 years"), format = "%Y")
lines(son.dates, m.index, col = "forestgreen", lwd = 2)
text(son.dates[25], s.index[25], "2006", pos = 3, cex = 0.9, col = "firebrick3")
text(son.dates[16], s.index[16], "1997", pos = 3, cex = 0.9, col = "firebrick3")
text(son.dates[13], s.index[13], "1994", pos = 3, cex = 0.9, col = "firebrick3")
text(son.dates[34], m.index[34], "2015", pos = 3, cex = 0.9, col = "forestgreen")
text(son.dates[6], m.index[6], "1987", pos = 3, cex = 0.9, col = "forestgreen")
text(son.dates[1], m.index[1], "1982", pos = 3, cex = 0.9, col = "forestgreen")
title("pIOD Indices (SON)", adj = 0)
abline(h=c(1.5, 1.25), lty = 2, lwd =1.5, col = c("firebrick3", "forestgreen"))
legend("topright",
       legend = c("S-Index", "M-Index"),
       lty    = 1,                 # line type
       lwd    = 2,                 # line width
       col    = c("firebrick3", "forestgreen"),
       #bty    = "n",               # no box; remove if you want a box
       inset  = 0.01)
dev.off()

#PCA plots
png(filename = "PCA_new.png", width = 2000, height = 2000, res = 300)
par(mar = c(5, 5, 3, 3))
plot(pc.std.IOD[,1], pc.std.IOD[,2], pch = 16, xlim = c(-2,4), ylim = c(-3, 5),
     xlab = "PC1", ylab = "PC2", col = "darkblue")
points(pc.std.IOD[c(13,16,25),1], pc.std.IOD[c(13,16,25),2], pch = 16, col = "firebrick3")
text(pc.std.IOD[c(13,16,25),1], pc.std.IOD[c(13,16,25),2], 
     c("1994", "1997", "2006"), pos = 4, cex = 0.9, col = "black")
points(pc.std.IOD[c(1,6,34),1], pc.std.IOD[c(1,6,34),2], pch = 16, col = "forestgreen")
text(pc.std.IOD[c(1,6,34),1], pc.std.IOD[c(1,6,34),2]+c(0.15, -0.1, 0), 
     c("1982", "1987", "2015"), pos = c(4,4,4), cex = 0.9, col = "black")
abline(h = 0.5, lty = 2)
abline(v = 1, lty = 2)
dev.off()


## strong years 1994, 1997, 2006 (index: 13, 16, 25)
## moderate years 1982, 1987, 2015 (index: 1, 6, 34)
## along with spatial mean, linear coeffs


#TODO: set-up preferred spatial colors
#match SST anomalies (deg C) from Cai et al 2021. n = 12, with most colors between (-1, 1)
#color match may be tricky since the lower and upper bounds may not be included
rybcol.12 <- rev(divergingx_hcl(n = 12, palette = "RdYlBu"))
rybcol.48 <- rev(divergingx_hcl(n = 48, palette = "RdYlBu"))

brks <- seq(eof.range2[1], eof.range2[2], length.out = length(roma_col) + 1)


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

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,1]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

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

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.diff[,,16]), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")#, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
world(add=TRUE)

#1993 SON anom detrended
k <- 13
temp <- array(sst.anom.avg, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.resid <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
brks <- seq(resid.range[1], resid.range[2], length.out = length(rybcol.48) + 1)
#dumb color attempt
cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))(12)
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
           col = rev(cols), breaks = breaks.dumb, zlim = resid.range, 
           xlab = "Lon", ylab = "Lat", #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
           axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
rect(40, -5, 100, 5, border = "black")
world(add=TRUE)

#just for iod domain
k <- 13
pIOD.resid <- aperm(sst.anom.pIOD, c(3, 2, 1)) 

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = pIOD.resid[,,k]), 
           col = rev(cols), breaks = breaks.dumb, zlim = resid.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
           #axes = FALSE)
world(add=TRUE)



##------TEST SECTION-----##
#TODO: delete when done
##test anomalies
A <- sst.array #as [lon, lat, time]

A.new <- aperm(A, perm = c(3, 2, 1))  # reorder data as [time, lat, lon]

#get dimensions from reordered array
nx <- dim(A.new)[3]
ny <- dim(A.new)[2] 
nt <- dim(A.new)[1]


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

#de-trend
#TODO: test linear de-trend and plot slope

#TODO: figure out how to work with irregular NAs in time/space (some years have different spatial NAs)
test.na <- which(!is.na(var), arr.ind = TRUE) 
true.index <- unique(test.na[,2])

A.resid <- matrix(NA, nrow = nt, ncol = ny*nx)
A.coef <- matrix(NA, nrow = 2, ncol = ny*nx)

var.x <- seq_len(nt) #time "covariates"
var.y <- var[ ,true.index] #response variable (sst anoms)

var.fit <- lm(var.y ~ var.x)

A.resid[,true.index] <- var.fit$residuals
A.coef[,true.index] <- var.fit$coefficients



#test viz (remove or move down when done):
#array reformatting code
temp <- array(var, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.diff <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]


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


# EOF test code
kmode <- 2
#sst.anom <- array(sst.anom.avg, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
sst.anom <- sst.anom.pIOD

#sst.anom as [time, lat, lon]
nt <- dim(sst.anom)[1]
ny <- dim(sst.anom)[2]
nx <- dim(sst.anom)[3]

#TODO: look at just importing as  [nt, ny*nx]
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

#TODO: double check negatives and splitting pc1 and pc2 for standarization
PC1.std <- scale(PC.temp[,1], center = TRUE, scale = TRUE)
PC2.std <- scale(PC.temp[,2], center = TRUE, scale = TRUE)
PC.std <- scale(PC.temp, center = TRUE, scale = TRUE)

#add in index plots 

#finalize the spatial eof (pca)
np <- nx*ny
V_eof <- matrix(NA, nrow = np, ncol = kmode)

#add in lsm sea data
V_eof[keep, ] <- EOF.temp