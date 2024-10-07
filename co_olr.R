## OLR code 

#TODO: move over working code from CO_OLR.rmd, or just start fresh with the aggregated data (olr_full.rda) 
## NOTE: we're starting fresh with aggregated data

#libraries
suppressMessages(library(lubridate))

suppressMessages(library(fields))

suppressMessages(library(scales))

#maping libraries

#suppressMessages(library(sp))
#suppressMessages(library(raster))
#suppressMessages(library(RColorBrewer))
#suppressMessages(library(terra))

# data
setwd("~/CO_AUS")
load("AUS_CO-main/olr_full.rda")


# additional data cleaning

## fix week 53 data 

temp_week53 <- olr_full_df[olr_full_df$Week == 53, ]
week53 <- which(olr_full_df$Week == 53, )


#loop through
for (k in week53) {
  olr_full_df[k-1, 1] <- (olr_full_df[k, 1] + olr_full_df[k-1, 1])/2
  olr_full_df[k-1, 2] <- (olr_full_df[k, 2] + olr_full_df[k-1, 2])/2
}

olr_full_df <- olr_full_df[-week53[3],]
olr_full_df <- olr_full_df[-week53[2],]
olr_full_df <- olr_full_df[-week53[1],]

rm(temp_week53, week53, k)

olr_full_df[olr_full_df$Week == 53, ]

# anomaly data

## weekly means

mean_df <- as.data.frame(matrix(NA, ncol=3))
colnames(mean_df) <- c("OLR_NE", "OLR_SE", "Week")

for (i in 1:52) {
  temp_week_df <- olr_full_df[olr_full_df$Week == i, ]
  temp_meanNE <- mean(temp_week_df$OLR_NE)
  temp_meanSE <- mean(temp_week_df$OLR_SE)
  new_row <- c(temp_meanNE, temp_meanSE, i)
  mean_df <- rbind(mean_df, new_row)
}

mean_df <- mean_df[-1, ]
row.names(mean_df) <- NULL

rm(i, temp_week_df, new_row, temp_meanNE, temp_meanSE)

## compute anoms

olr_anom_df <- as.data.frame(matrix(NA, ncol = 4))
colnames(olr_anom_df) <- c("NE_anom", "SE_anom", "Date", "Week")

length(olr_full_df$Date)

for (j in 1:length(olr_full_df$Date)) {
  temp_week <- olr_full_df[j, 4]
  NE_mean <- mean_df[temp_week, 1]
  SE_mean <- mean_df[temp_week, 2]
  
  NE_anom <- olr_full_df[j, 1] - NE_mean
  SE_anom <- olr_full_df[j, 2] - SE_mean
  
  temp_date <- olr_full_df[j, 3]
  
  temp_row <- c(NE_anom, SE_anom,  as_date(temp_date), temp_week)
  olr_anom_df <- rbind(olr_anom_df, temp_row)
}

olr_anom_df <- olr_anom_df[-1, ]
olr_anom_df$Date <- as_date(olr_anom_df$Date)
row.names(olr_anom_df) <- NULL

rm(j, temp_week, NE_mean, SE_mean, temp_date, NE_anom, SE_anom, temp_row)


setwd("~/CO_AUS/Aus_CO-main")
save(olr_anom_df, file = "olr_anom.rda")

#test section (delete when done)



# plots

## basic plots

png(filename = "OLR_mean.png", width = 2200, height = 1200, res = 300)
  plot(mean_df$Week, mean_df$OLR_NE, type = "l", col= "red", 
     ylim = c(200,300), main = "Weekly OLR Mean (2000 to 2019)", xlab = "Week", ylab = "OLR [w/m^2]")
  lines(mean_df$Week, mean_df$OLR_SE, col = "blue")
  legend("bottomright", inset = c(0, 0), legend = c("NE Aus", "SE Aus"), 
       col = c("red", "blue"), lty = c(1,1))
dev.off()


## TS plots (from CO_base_vis.rmd)

### setup
olr_anoms <- olr_anom_df

top.color <- "orange"
bottom.color <- "cyan"

names(olr_anoms) <- c( "NE Aus", "SE Aus", "time", "week")

y2.lab.vals <- c("Anomaly OLR [w/m^2]",
                 "Anomaly OLR [w/m^2]")

olr <- olr_anoms

olr_names <- names(olr)

#normalizing for both values
for (i in 1:2) {
  olr[,i] <- (olr[,i] - mean(olr[,i]))/ sd(olr[,i])
}

first_date <- min(olr$time)
last_date <- max(olr$time)

time <- olr$time
time_range <- c(first_date, last_date)
x_ticks <- seq(year(time_range[1]), year(time_range[2]), by = 1) #can also be done with lubridate floor
x_ticks <- ymd(paste0(x_ticks, "01", "01"))

### plot (wrap with png when done)

png(filename = "OLR_ts.png", width = 2500, height = 2000, res = 300)

par(mfrow = c(2, 1))
par(mar = c(0,2,0,0))
par(oma = c(1.5,2,0,0))
par(mgp = c(2.75,1,0))

for (i in 1:2){
  this.olr <- (olr[,i])
  temp.time <- as.Date(olr$time)
  #temp.time <- gsub("-", "", temp.time)
  
  y.tick.max <- max(round(range(this.olr)))
  # y.tick.max <- 4
  y.tick.max.lab <- y.tick.max / 2
  y.ticks <- seq(-y.tick.max, y.tick.max, by = 1)
  y.tick.labs <- c(-y.tick.max.lab, 0, y.tick.max.lab)
  
  plot(temp.time, this.olr, type = "l", col = "black", lwd = 2,
       xaxt = "n", xlab = "",
       yaxt = "n", ylab = y2.lab.vals[i], col.lab = "black",
       xlim = c(as.Date(time_range[1]) + months(7), as.Date(time_range[2]) - months(7)),
       ylim = range(y.ticks), bty = "n", cex.lab = 1.25,  xpd = NA)
  
  if (i == 2){
    text(x = x_ticks + months(6),
         y = range(y.ticks)[1]-0.7,
         labels = year(x_ticks),
         cex = 1.3, col = "black", xpd = NA)
  }
  
  
  if (i == 1){ abline(v = x_ticks[1:(length(x_ticks))],
                      lty = 2, col = "grey", lwd = 2, xpd = NA) }
  
  axis(side = 2, at = y.tick.labs, labels = y.tick.labs, col = NA, 
       col.ticks = "black", col.axis = "black")
  
  abline(h = 0, lty = 1, col = "grey", lwd = 1)
  
  over <- this.olr >= 0
  this.olr.over <- this.olr
  this.olr.over[!over] <- 0
  this.olr.under <- this.olr
  this.olr.under[over] <- 0
  
  envelopePlot(x1 = temp.time,
               y1 = this.olr.over,
               x2 = temp.time,
               y2 = rep(0, length(time)),
               col = alpha(top.color, 0.7),
               lineCol = NA)
  envelopePlot(x1 = temp.time,
               y1 = this.olr.under,
               x2 = temp.time,
               y2 = rep(0, length(time)),
               col = alpha(bottom.color, 0.7),
               lineCol = NA)
  
  
  #td <- y.tick.labs[3]
  #legend(x = c(ymd("1999-05-01"), ymd("2001-10-01")),
  #       y = c(2.0*td, 3.25*td),
         # legend = toupper(pred.names[i]),
  #       legend = olr_names[i],
  #       box.col = "white", bg = "white",
         # adj = c(0., 0.6),
  #      xpd = NA, text.col = "grey40", cex = 1.5)
  
}
dev.off()
