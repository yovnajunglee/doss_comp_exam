### Process data
rm(list=ls())
library(tidyverse)
library(dplyr)
# Load data
load("~/Desktop/Research Comps/Data/2017/2017_max_temp_data.RData") 
load("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tempmax_grid_60km_2017.RData") 
source("~/Code/GaussianTempModel.R")
source("~/Desktop/Research Comps/Code/Author/Personal/GaussianTempModel.R")

vecz1 = vecz
time_out = 266:275
vecz = vecz1[-time_out,]
J= nrow(vecz) ;I = ncol(vecz) ### 303 por 61
z = vecz
Inew = I


lat<- info.sites4$lat ; long<- info.sites4$long ; 
x1<-info.sites4$type
Plat<- (lat - median(lat))/diff(quantile(lat, c(0.25,0.75))) 
Plong<- (long - median(long))/diff(quantile(long, c(0.25,0.75))) 
Ptemp = (t(dt.temp4[,-time_out])-median(t(dt.temp4[,-time_out])))/diff(quantile(dt.temp4[,-time_out], c(0.25,0.75))) 


# Out-of-sample time-points for forecasts at observed locations 
Ptemp.p = (t(dt.temp4[,time_out])-median(t(dt.temp4[,time_out])))/diff(quantile(dt.temp4[,time_out], c(0.25,0.75))) 


# Out-of-sample locations for spatial predictions
oos_data = tempmax_grid_2017 %>% filter(month >=3, month <=11)
info.sites.out = oos_data %>% dplyr::select(long, lat) %>% distinct() %>% mutate(code = 1:n())
oos_data <- left_join(oos_data, info.sites.out, by=c("long","lat"))
oos_data <- oos_data %>% arrange(code, month, day)
nout = oos_data %>% dplyr::select(long, lat) %>% distinct() %>% nrow()
lat.out<- info.sites.out$lat ; long.out<- info.sites.out$long; 
Plat.out<- (lat.out - median(lat.out))/diff(quantile(lat.out, c(0.25,0.75))) 
Plong.out<- (long.out - median(long.out))/diff(quantile(long.out, c(0.25,0.75)))
dt.tempout = matrix(oos_data$tempmax, nrow = nout, byrow = T)
Ptemp.out = (t(dt.tempout[,-time_out])-median(t(dt.tempout[,-time_out])))/diff(quantile(dt.tempout[,-time_out], c(0.25,0.75))) 

## coordinates
coords<- cbind(lat,long)

### distance matrix
distS<- distf(coords)
median(distS)
max(distS)

### distance matrix
coords.out<- cbind(lat.out,long.out)
distS.out <-  distf(coords.out)
distS.tot<- as.matrix(distf(rbind(coords,coords.out)))

#==== Fit Models


# 1) Fit Gaussian model with nugget effect
source("~/Desktop/Research Comps/Code/Author/Personal/run_gaussian_temp.R")
gaussian_grid <- run_gaussian_temp(M=10000)

source("~/Desktop/Research Comps/Code/Author/Personal/run_covdynglg_temp.R")
covdynglg_grid <- run_covdynglg_temp(M=20000)

source("~/Desktop/Research Comps/Code/Author/Personal/run_covdynglg_matern_temp.R")
covdynglg_matern_grid_2<- run_covdynglg_matern(M=20e3, nu=.5)


