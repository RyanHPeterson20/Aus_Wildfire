
#R script for finding SST anomalies and EOF/PCAs

#libraries
#.nc files
suppressMessages(library(ncdf4))
suppressMessages(library(terra))

# date/data mgmt
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
load("sst_extend.rda")


#setup 
lsm.IOD.array <- array(lsm.IOD, dim = dim(sst.OISST.new))

sst.OISST.masked <- ifelse(lsm.IOD.array == 1, sst.OISST.new, NA)

#updated (for up to 2019)
lsm.IOD.array2 <- array(lsm.IOD, dim = dim(sst.OISST.new2))

sst.OISST.masked2 <- ifelse(lsm.IOD.array2 == 1, sst.OISST.new2, NA)

#dim check
dim(sst.OISST.masked2)
dim(sst.godas2)

#TODO: check for NA's across the entire period (that is, find a uniform LSM)
#check on SODA
##which(is.na(sst.interp.SODA))
#Check on GODAS

#internal functions


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


##----Main----##


#sst anomalies SON, for 
##OISST
sst.anom.oisst <- sst.anoms(sst.OISST.masked)
#GODAS
sst.anom.godas <- sst.anoms(sst.array)
#SODA
sst.anom.soda <- sst.anoms(sst.interp.SODA)

#get SON average of SST anomalies
sst.anom.avg <- (sst.anom.oisst$anom + sst.anom.godas$anom + sst.anom.soda$anom)/3

#TODO: add back in pca for "full" region
nt <- 34
nx <- length(lon.values)
ny <- length(lat.values)

sst.anom <- array(sst.anom.avg, dim = c(nt, ny, nx)) 
pca.base <- sst.eof(sst.anom, kmode = 2)

eof.base <- pca.base$EOF


#select for reduced region pIOD (redo from earlier, redundant but oh-well)
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

V_eof <- pca.pIOD$EOF
pc.std.IOD <- scale(pca.pIOD$PC, center = TRUE, scale = TRUE)

pc.std.IOD[,1] <- -pc.std.IOD[,1] 

PC1 <- pc.std.IOD[,1]
PC2 <- pc.std.IOD[,2]

s.index <- (pc.std.IOD[,1] + pc.std.IOD[,2])/sqrt(2)
m.index <- (pc.std.IOD[,1] - pc.std.IOD[,2])/sqrt(2)

#TODO: save these further down
plot(1:34, s.index, type = "l", col = "firebrick", ylim=c(-2,5))
lines(1:34, m.index, col = "darkgreen")
abline(h=c(1.5, 1.25), lty = 2, col = c("firebrick", "darkgreen"))

#PC-ts for PC1 and PC2
plot(1:34, pc.std.IOD[,1], type = "l")
plot(1:34, pc.std.IOD[,2], type = "l")

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


#--sst regression to here: (fig 1a and 1b)
#project SST anoms onto each EOF pattern
eof.mask <- is.finite(eof.base[,1]) #currently using eof1

v.eof1 <- eof.base[eof.mask,1] #eof1 vector (may or may NOT use this)
#TODO: potentially average here? ("for the 1982-2015 period average" - cai et al)
x.anoms <- sst.anom.avg[,eof.mask] #sst anom data matrix

spatial.proj <- matrix(NA, ncol = 2, nrow = 3600)
proj1 <- t(x.anoms) %*% pc.std.IOD[,1]
proj2 <- t(x.anoms) %*% pc.std.IOD[,2]

#proj1 <- t(x.anoms) %*% pca.pIOD$PC[,1]
#proj2 <- t(x.anoms) %*% pca.pIOD$PC[,2]

spatial.proj[eof.mask, ] <- cbind(proj1, proj2)


temp <- array(t(spatial.proj), dim = c(2, ny, nx))  # t(var) is [nt, ny*nx]
#v.eof1 <- temp[1,,]
#v.eof2 <- temp[2,,]
proj.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]

proj.spatial <- proj.spatial/34


#TODO: finalize WIP for figs 1e and 1f
#work on this later, I might be onto something
EOF1.full <- proj.spatial[,,1]
EOF2.full <- proj.spatial[,,2]

test.1 <- PC1 %*% t(proj1) 
test.2 <- PC2 %*% t(proj2)

spatial.test1 <- matrix(NA, ncol = 3600, nrow = 34)
spatial.test2 <- matrix(NA, ncol = 3600, nrow = 34)

spatial.test1[ ,eof.mask] <- test.1
spatial.test2[ ,eof.mask] <- test.2


#test location [6,1]
proj.spatial[6,1,1] +proj.spatial[6,1,2]
s.eof[6,1]

imagePlot(list(x = lon.values, y = rev(lat.values), z = s.eof), 
      xlab = "Lon", ylab = "Lat")
world(add=TRUE)

#setup color and range
eof.absmax <- max(abs(proj.spatial), na.rm = TRUE)
eof.range <- c(-eof.absmax, eof.absmax) 

#color setup
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
brks1 <- c(eof.range[1], seq(-0.5, 0.5, 0.1), eof.range[2])
brks <- seq(eof.range[1], eof.range[2], length.out = length(cols.ryb) + 1)
leg.breaks <- seq(-0.6, 0.6, 0.1)


setwd("~/CO_AUS/Aus_CO-main/pIOD/Figures")

#figs 1a and 1b reproduced here.
png(filename = "eof1_1a.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = proj.spatial[,,1]), 
      col = cols.ryb, breaks = brks1, zlim = eof.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-0.6, 0.6), legend.only = TRUE, col = cols.ryb,
           breaks = leg.breaks, horizontal = TRUE, 
           axis.args = list(at = seq(-0.5, 0.5, 0.1)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
dev.off()

png(filename = "eof2_1b.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = proj.spatial[,,2]), 
      col = cols.ryb, breaks = brks1, zlim = eof.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-0.6, 0.6), legend.only = TRUE, col = cols.ryb,
           breaks = leg.breaks, horizontal = TRUE, 
           axis.args = list(at = seq(-0.5, 0.5, 0.1)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
dev.off()


#transparent saves
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
temp <- array(t(V_eof), dim = c(2, 10, 60))  # t(var) is [nt, ny*nx]
#v.eof1 <- -temp[1,,]
#v.eof2 <- temp[2,,]
eof.spatial <- aperm(temp, c(3, 2, 1))  # [nx, ny, nt]
#eof.spatial <- aperm(Q_net_1, c(2, 1, 3))

#setup color and range
#setup color and range
eof.absmax <- max(abs(eof.spatial), na.rm = TRUE)
eof.range <- c(-eof.absmax, eof.absmax) 

#color setup
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
brks <- seq(eof.range[1], eof.range[2], length.out = length(cols.ryb) + 1)

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = -eof.spatial[,,1]), 
           col = cols.ryb, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = eof.spatial[,,2]), 
           col = cols.ryb, breaks = brks, zlim = eof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

#naive s.eof, m.eof
eof1.std <- scale(-eof.spatial[,,1])
eof2.std <- scale(eof.spatial[,,2])
s.eof <- (eof1.std + eof2.std)/sqrt(2)
m.eof <- (eof1.std - eof2.std)/sqrt(2)

pIODeof.absmax <- max(abs(s.eof), abs(m.eof), na.rm = TRUE)
pIODeof.range <- c(-pIODeof.absmax, pIODeof.absmax) 
brks.new <- seq(pIODeof.range[1], pIODeof.range[2], length.out = length(cols.ryb) + 1)

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = s.eof), 
           col = cols.ryb, breaks = brks.new, zlim = pIODeof.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = m.eof), 
           col = cols.ryb, breaks = brks.new, zlim = pIODeof.range, 
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

##SST anomalies (detrended)
temp <- array(sst.anom.avg, dim = c(nt, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.resid <- aperm(temp, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

#1994 SON anom detrended
k <- 13
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 


setwd("~/CO_AUS/Aus_CO-main/pIOD/Figures")

#extended fig1a
png(filename = "sst1994anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
           col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
           xlab = "Lon", ylab = "Lat",
           axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("1994: SST Anomalies", adj = 0)
dev.off()


#1997 SON anom detrended
k <- 16
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 

#extended fig1b
png(filename = "sst1997anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
      col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("1997: SST Anomalies", adj = 0)
dev.off()

#2006 SON anom detrended
k <- 25
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 

#extended fig1c
png(filename = "sst2006anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
      col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("2006: SST Anomalies", adj = 0)
dev.off()


#1982 SON anom detrended
k <- 1
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 

#extended fig2a
png(filename = "sst1982anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
      col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("1982: SST Anomalies", adj = 0)
dev.off()


#1987 SON anom detrended
k <- 6
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 

#extended fig2b
png(filename = "sst1987anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
      col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("1987: SST Anomalies", adj = 0)
dev.off()


#2015 SON anom detrended
k <- 34
resid.absmax <- max(abs(spat.resid[,,k]), na.rm = TRUE)
resid.range <- c(-resid.absmax, resid.absmax) 
#default breaks
cols.ryb <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(12))
breaks.dumb <- c(resid.range[1], seq(-1,1, 0.2), resid.range[2])
leg.brks2 <- seq(-1.2, 1.2, 0.2) 

#extended fig2c
png(filename = "sst2015anom.png", width = 2400, height = 1600, res = 250)
par(mar = c(7, 5, 2, 5), mgp = c(2.25, 1, 0))
image(list(x = lon.values, y = rev(lat.values), z = spat.resid[,,k]), 
      col = cols.ryb, breaks = breaks.dumb, zlim = resid.range, 
      xlab = "Lon", ylab = "Lat",
      axes = FALSE)
axis(1, at = c(60, 90, 120))  
axis(2)                      
box()
world(add=TRUE)
image.plot(zlim = c(-1.2, 1.2), legend.only = TRUE, col = cols.ryb,
           breaks = leg.brks2, horizontal = TRUE, 
           axis.args = list(at = seq(-1, 1, 0.2)),
           smallplot = c(0.1, 0.9, 0.08, 0.10))
title("2015: SST Anomalies", adj = 0)
dev.off()


#Get spatial mean and coefs for 1982 to 2015
#godas mean
spat.mean.godas <- array(sst.anom.godas$mean, dim = c(ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.mean.oisst <- array(sst.anom.oisst$mean, dim = c(ny, nx))

image.plot(list(x = lon.values, y = rev(lat.values), z = t(spat.mean.oisst)), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)


#Get spatial mean and coefs for 1982 to 2019
sst.anom.godas2019 <- sst.anoms(sst.godas2)
sst.anom.oisst2019 <- sst.anoms(sst.OISST.masked2)

spat.mean.godas2019 <- array(sst.anom.godas2019$mean, dim = c(ny, nx))
spat.mean.oisst2019 <- array(sst.anom.oisst2019$mean, dim = c(ny, nx))

image.plot(list(x = lon.values, y = rev(lat.values), z = t(spat.mean.oisst2019)), 
           col = tim.colors(256), 
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

#get difference between the two
oisst.diff <- (spat.mean.oisst2019 - spat.mean.oisst)
godas.diff <- (spat.mean.godas2019 - spat.mean.godas)

diff.max <- max(abs(oisst.diff), na.rm = TRUE)+0.02
diff.range <- c(-diff.max , diff.max)
#TODO: get red-blue diverge colors
cols.rb <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(48))

image.plot(list(x = lon.values, y = rev(lat.values), z = t(oisst.diff)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = t(godas.diff)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

#Get spatial mean and coefs for 1982 to 2018
sst.anom.godas2018 <- sst.anoms(sst.godas2[,,1:444])
sst.anom.oisst2018 <- sst.anoms(sst.OISST.masked2[,,1:444])

spat.mean.godas2018 <- array(sst.anom.godas2018$mean, dim = c(ny, nx))
spat.mean.oisst2018 <- array(sst.anom.oisst2018$mean, dim = c(ny, nx))
#get difference between the two
oisst.diff2018 <- (spat.mean.oisst2018 - spat.mean.oisst)
godas.diff2018 <- (spat.mean.godas2018 - spat.mean.godas)

image.plot(list(x = lon.values, y = rev(lat.values), z = t(oisst.diff2018)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = t(godas.diff2018)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)


#2019-2018
oisst.difflate <- (spat.mean.oisst2019 - spat.mean.oisst2018)
godas.difflate <- (spat.mean.godas2019 - spat.mean.godas2018)

image.plot(list(x = lon.values, y = rev(lat.values), z = t(oisst.difflate)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = t(godas.difflate)), 
           col = cols.rb, zlim = diff.range,
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)


#Coeffs from 1982-2015 and 2019

## oisst/godas coefs from 1982-2015
#coef.max <- max(abs(sst.anom.godas$coef), na.rm = TRUE)
coef.max <- max(abs(sst.anom.oisst$coef), na.rm = TRUE)

temp1 <- array(sst.anom.godas$coef, dim = c(2, ny, nx))  # from [nt, ny*nx] to [nt, ny, nx]
spat.coefs.godas <- aperm(temp1, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

temp2 <- array(sst.anom.oisst$coef, dim = c(2, ny, nx)) 
spat.coefs.oisst <- aperm(temp2, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

range(spat.coefs[,,2], na.rm = TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,1]), 
           col = cols.rb, zlim = c(-coef.max, coef.max),
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,2]), 
           col = cols.rb, zlim = c(-0.075, 0.075),
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

#coeffs from 1982-2019
coef.max2019 <- max(abs(sst.anom.godas2019$coef), na.rm = TRUE)

temp1 <- array(sst.anom.godas2019$coef, dim = c(2, ny, nx)) 
spat.coefs.godas2019 <- aperm(temp1, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]

temp2 <- array(sst.anom.oisst2019$coef, dim = c(2, ny, nx)) 
spat.coefs.oisst2019 <- aperm(temp2, c(3, 2, 1))  # from [nt, ny, nx] to [nx, ny, nt]


image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,1]), 
           col = cols.rb, zlim = c(-coef.max2019, coef.max2019),
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)

image.plot(list(x = lon.values, y = rev(lat.values), z = spat.coefs[,,2]), 
           col = cols.rb, zlim = c(-0.075, 0.075),
           xlab = "Lon", ylab = "Lat")
world(add=TRUE)


#coeff differences






##-------- TEST SECTION -------##


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
dev.off()


#just for iod domain
k <- 13
pIOD.resid <- aperm(sst.anom.pIOD, c(3, 2, 1)) 

image.plot(list(x = lon.values.IOD, y = rev(lat.values.IOD), z = pIOD.resid[,,k]), 
           col = rev(cols), breaks = breaks.dumb, zlim = resid.range, 
           xlab = "Lon", ylab = "Lat") #, main = paste0("OISST v2: ", format(as.Date(times[25+1]), "%B %Y")))
#axes = FALSE)
world(add=TRUE)


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