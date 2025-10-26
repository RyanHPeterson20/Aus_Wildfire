#data preparation/cleaning for Australian CO Modeling (wildfire/climate modes)


#TODO: bring over existing data cleaning and merge it into a single file/script



#library
suppressMessages(library(ncdf4)) #.nc files
suppressMessages( library(lubridate))

#load in functions (as needed)

#raw data import (from .csv)

#import from .nc
#TODO: update the source .nc files
setwd("~/CO_AUS/Aus_CO-main")
#DMI
dmi.nc <- nc_open("~/CO_AUS/Aus_CO-main/pIOD/dmi(2).nc")
#WTIO
wtio.nc <- nc_open("~/CO_AUS/Aus_CO-main/pIOD/wtio(1).nc")
#ETIO
etio.nc <- nc_open("~/CO_AUS/Aus_CO-main/pIOD/setio(2).nc")




## --- Main --- ##



## IOD: (DMI, WTIO, ETIO)

#weeks (time)
time.dmi <- ncvar_get(dmi.nc, "WEDCEN2") # days since 1800-01-01 00:00:0.0
time.dmi <- as_datetime("1900-01-01T00:00:00") + days(time.dmi)
time.wtio <- ncvar_get(wtio.nc, "WEDCEN2") # days since 1800-01-01 00:00:0.0
time.wtio <- as_datetime("1900-01-01T00:00:00") + days(time.wtio)
time.etio <- ncvar_get(etio.nc, "WEDCEN2") # days since 1800-01-01 00:00:0.0
time.etio <- as_datetime("1900-01-01T00:00:00") + days(time.etio)

#weekly data
sst.dmi <- ncvar_get(dmi.nc, "DMI")
sst.wtio <- ncvar_get(wtio.nc, "WTIO")
sst.etio <- ncvar_get(etio.nc, "SETIO")

#get df
dmi.new.df <- data.frame(sst = sst.dmi, date = time.dmi, week = epiweek(time.dmi))
wtio.df <- data.frame(sst = sst.wtio, date = time.wtio, week = epiweek(time.wtio))
etio.df <- data.frame(sst = sst.etio, date = time.etio, week = epiweek(time.etio))

#from 1982-01-06 to 2020-04-01
time.domain <- which(dmi.new.df$date >= time.dmi[10] & dmi.new.df$date <= time.dmi[2005])

dmi.new1.df <- dmi.new.df[time.domain,]
dmi.new1.df$year <- year(dmi.new1.df$date)
dmi.new1.df$month <- month(dmi.new1.df$date)

wtio1.df <- wtio.df[time.domain,]
wtio1.df$year <- year(wtio1.df$date)
wtio1.df$month <- month(wtio1.df$date)

etio1.df <- etio.df[time.domain,]
etio1.df$year <- year(etio1.df$date)
etio1.df$month <- month(etio1.df$date)



