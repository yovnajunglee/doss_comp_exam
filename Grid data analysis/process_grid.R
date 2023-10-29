library(raster)
#install.packages("ncdf4")
library(ncdf4)
library(ggplot2)
library(tidyverse)


## Process grid data
# ** Function to create daily grid
create_daily_grid <- function(long, lat,  temp_array, n){
  I = temp_array[,,1]*0 + 1
  dat <- data.frame(day = 1, long = long, lat = lat, tempmax = c(temp_array[,,1]))
  for(i in 2:n){
    dat <- rbind(dat,data.frame(day = i, long = long, lat = lat, tempmax = c(temp_array[,,1])) )
  }
  return(na.omit(dat))
}


## Jan 2017
r1 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170101-20170131.nc")
lat1 = ncvar_get(r1, "latitude" )
lon1 = ncvar_get(r1, "longitude")
n1 = dim(ncvar_get(r1, "time" ))
temp_array1 <- ncvar_get(r1, "tasmax")
jan_2017 <- create_daily_grid(c(lon1), c(lat1), temp_array1, n1 )




## Feb 2017
r2 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170201-20170228.nc")
lat2 = ncvar_get(r2, "latitude" )
lon2 = ncvar_get(r2, "longitude")
n2 = dim(ncvar_get(r2, "time" ))
temp_array2<- ncvar_get(r2, "tasmax")
feb_2017 <- create_daily_grid(c(lon2), c(lat2), temp_array2, n2 )


## March 2017
r3 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170301-20170331.nc")
lat3 = ncvar_get(r3, "latitude" )
lon3 = ncvar_get(r3, "longitude")
n3 = dim(ncvar_get(r3, "time" ))
temp_array3<- ncvar_get(r3, "tasmax")
march_2017 <- create_daily_grid(c(lon3), c(lat3), temp_array3, n3 )



## April 2017
r4 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170401-20170430.nc")
lat4 = ncvar_get(r4, "latitude" )
lon4 = ncvar_get(r4, "longitude")
n4 = dim(ncvar_get(r4, "time" ))
temp_array4<- ncvar_get(r4, "tasmax")
april_2017 <- create_daily_grid(c(lon4), c(lat4), temp_array4, n4 )

           
## May 2017
r5 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170501-20170531.nc")
lat5 = ncvar_get(r5, "latitude" )
lon5 = ncvar_get(r5, "longitude")
n5 = dim(ncvar_get(r5, "time" ))
temp_array5<- ncvar_get(r5, "tasmax")
may_2017 <- create_daily_grid(c(lon5), c(lat5), temp_array5, n5 )

## June 2017
r6 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170601-20170630.nc")
lat6 = ncvar_get(r6, "latitude" )
lon6 = ncvar_get(r6, "longitude")
n6 = dim(ncvar_get(r6, "time" ))
temp_array6<- ncvar_get(r6, "tasmax")
june_2017 <- create_daily_grid(c(lon6), c(lat6), temp_array6, n6 )


## July 2017
r7 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170701-20170731.nc")
lat7 = ncvar_get(r7, "latitude" )
lon7 = ncvar_get(r7, "longitude")
n7 = dim(ncvar_get(r7, "time" ))
temp_array7<- ncvar_get(r7, "tasmax")
july_2017 <- create_daily_grid(c(lon7), c(lat7), temp_array7, n7 )


## Aug 2017
r8 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170801-20170831.nc")
lat8 = ncvar_get(r8, "latitude" )
lon8 = ncvar_get(r8, "longitude")
n8 = dim(ncvar_get(r8, "time" ))
temp_array8<- ncvar_get(r8, "tasmax")
aug_2017 <- create_daily_grid(c(lon8), c(lat8), temp_array8, n8 )

## Sept 2017
r9 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20170901-20170930.nc")
lat9 = ncvar_get(r9, "latitude" )
lon9 = ncvar_get(r9, "longitude")
n9 = dim(ncvar_get(r9, "time" ))
temp_array9<- ncvar_get(r9, "tasmax")
sept_2017 <- create_daily_grid(c(lon9), c(lat9), temp_array9, n9 )


## Oct 2017
r10 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20171001-20171031.nc")
lat10 = ncvar_get(r10, "latitude" )
lon10 = ncvar_get(r10, "longitude")
n10 = dim(ncvar_get(r10, "time" ))
temp_array10<- ncvar_get(r10, "tasmax")
oct_2017 <- create_daily_grid(c(lon10), c(lat10), temp_array10, n10 )

## Nov 2017
r11 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20171101-20171130.nc")
lat11 = ncvar_get(r11, "latitude" )
lon11 = ncvar_get(r11, "longitude")
n11 = dim(ncvar_get(r11, "time" ))
temp_array11<- ncvar_get(r11, "tasmax")
nov_2017 <- create_daily_grid(c(lon11), c(lat11), temp_array11, n11 )


# Dec 2017
r12 = nc_open("~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tasmax_hadukgrid_uk_60km_day_20171201-20171231.nc")
lat12 = ncvar_get(r12, "latitude" )
lon12 = ncvar_get(r12, "longitude")
n12 = dim(ncvar_get(r12, "time" ))
temp_array12<- ncvar_get(r12, "tasmax")
dec_2017 <- create_daily_grid(c(lon12), c(lat12), temp_array12, n12 )


nlocs =  78#nrow(feb_2017)/28





tempmax_grid_2017  = rbind(cbind(month =1, jan_2017),
                           cbind(month =2, feb_2017),
                           cbind(month =3, march_2017),
                           cbind(month =4, april_2017),
                           cbind(month =5, may_2017),
                           cbind(month =6, june_2017),
                           cbind(month =7, july_2017),
                           cbind(month =8, aug_2017),
                           cbind(month =9, sept_2017),
                           cbind(month =10, oct_2017),
                           cbind(month =11, nov_2017),
                           cbind(month =12, dec_2017))

save(tempmax_grid_2017, file="~/Desktop/Research Comps/Data/2017/UK_Grid_TempMax/tempmax_grid_60km_2017.RData")

worldmap = map_data('world')
ggplot() +
  geom_polygon(data = worldmap  %>% filter(region == "UK"),
               aes(x = long,
                   y = lat,
                   group = group), fill = "white", colour = 1) +
  #geom_point(data = tempmax_grid_2017 %>% filter(day==30, month==8), aes(x=long, y=lat, col = tempmax))+
  geom_point(data = info.sites4, aes(x = long, y =lat), col =2, shape=3)+
  geom_point(aes(y = 54.56930,
                 x = -1.220874)) + geom_text(data = info.sites4, aes(label = code, x = long, y =lat))

