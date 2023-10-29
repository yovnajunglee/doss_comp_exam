## UK Ozone data


#install.packages("openair")
library(openair)
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(worldmet)
library(ggcorrplot)
library(gstat)
library(spacetime)
library(randomForest)
library(missForest)


## Extract station metadata

info_sites <- importMeta(
  all = TRUE,
  year = 2017
)

## Extract air pollution (ozone) data 
## for  2017

sites_max_o3_2017 <- importAURN(
  year = 2017,
  site = info_sites$code,
  pollutant = "o3",
  data_type = "daily_max_8",
  meta = TRUE
)
names(sites_max_o3_2017)
save(sites_max_o3_2017, file = "sites_max_o3_2017.RData")

## Data cleaning
#sites_pollution_2017$day <- substring(sites_pollution_2017$date,1,10)
sites_max_o3_2017$month <- as.numeric(substring(sites_max_o3_2017$date,6,7))

## Remove duplicate rows
max_o3_2017 <- sites_max_o3_2017 %>%  distinct()

## No. of days
ndays = length(unique(max_o3_2017$date));ndays
site_in <- (max_o3_2017 %>% group_by(code) %>% summarise(n = n()) %>% filter(n == ndays) );site_in
#max_o3_2017 <- max_o3_2017 %>% filter(code %in% site_in$code)


## Extract meteorological data 
## for November 2017
meteo_2017 <- importAURN(
  year = 2017,
  site = info_sites$code,
  meta = TRUE
)
meteo_2017$day <- substring(meteo_2017$date,1,10)
save(meteo_2017, file = "meteo_2017.RData")
names(meteo_2017)

# Find average temperature per day
daily_mean_2017 <- meteo_2017 %>% dplyr::select(code, day, air_temp, ws, o3) %>%
  group_by(code, day) %>%
  summarise(mean_daily_temp = mean(air_temp), 
            mean_daily_ws = mean(ws))
max_o3_2017$day <- as.character(max_o3_2017$date)
all_data_2017 <- left_join(max_o3_2017, daily_mean_2017, by = c("code" = "code", 
                                                                "day" = "day")) 

## Only select stations with < 5% missing data 
## Extract data for March - Nov 
all_data_2017 <-  all_data_2017 %>% filter(code %in% summ_missing$code, month >= 3, month <= 11)

length(unique(all_data_2017$code))
length(unique(all_data_2017$site))
length(unique(all_data_2017$date))

## Preliminary EDA
all_data_2017 %>% ggplot(aes(x = date, y = o3)) + geom_line()
all_data_2017 %>% ggplot(aes(x = date, y = mean_daily_temp)) + geom_line()
all_data_2017 %>% ggplot(aes(x = date, y = mean_daily_ws)) + geom_line()

names(all_data_2017)

#save(all_data_2017, file = "all_data_2017.RData")

## Create data matrix

#load("all_data_2017.RData")
# Remove stations not in paper
not_inc <- c("MH", "DERR", "LN", "BEL2")
### Combined data (2017 and 2018)
#all_data <- rbind(cbind(all_data_2017, year =2017),cbind(all_data_2018, year = 2018)) 
all_data <- all_data_2017
all_data <- all_data %>% arrange(code, date)

## Find stations with < 5% missing data
ndays = length(unique(all_data$date))
code <- unique(all_data$code)
days <- data.frame(day = unique(all_data$day))
summ_missing <- all_data  %>% group_by(code) %>% 
  summarise(max_8h_o3 = sum(is.na(o3))*100/ndays,
            mean_daily_temp = sum(is.na(mean_daily_temp))*100/ndays,
            mean_daily_ws =  sum(is.na(mean_daily_ws))*100/ndays)
summ_missing <- summ_missing %>% 
  filter(max_8h_o3 < 5 , mean_daily_temp < 5 , mean_daily_ws < 5)
nstations = nrow(summ_missing)
nstations*ndays

## Filter out stations with > 5% missing data
all_data <-  all_data %>% filter(code %in% summ_missing$code)
site_out <- all_data %>% group_by(code) %>% 
  summarise(n = n()) %>% filter(n!=ndays)
all_data1 <- all_data %>% filter(!code%in%c(not_inc, site_out$code))
nstations1 = length(unique(all_data1$code))
code1 = unique(all_data1$code)

#which(code1 %in% c("BORN","NOTT"))

all_data1%>%gplot(aes(x=date,y=mean_daily_temp, col = code))+geom_line()

# when using 2017-2018 index is 6,37
nstations1*ndays


## Create data matrix

dt.temp4 <- (matrix(all_data1$mean_daily_temp, byrow = T,  nrow=nstations1))
dt.ws4 <- (matrix(all_data1$mean_daily_ws, byrow = T,  nrow=nstations1))
vecz <- t(matrix(all_data1$o3, byrow = T, nrow = nstations1))

## Impute data using a random forest algorithm
o3_mis <- missForest(vecz) 
vecz = o3_mis$ximp

temp_mis <- missForest(t(dt.temp4)) 
dt.temp4 = t(temp_mis$ximp)


ws_mis <- missForest(t(dt.ws4)) 
dt.ws4 = t(ws_mis$ximp)

info.sites4 <- all_data1 %>%
  dplyr::select(code, latitude, longitude) %>% distinct()
colnames(info.sites4)[2:3] <- c("lat", "long")


