library(sf)
library(stars)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw(base_size = 16))


## Create matrix of predictions
M=20e3
nstation = 78
predicted_gaussian = matrix(0, nrow = 78*265, ncol = 1)
predicted_covdynglg = matrix(0, nrow = 78*265, ncol = 1)
predicted_covdynglg_matern = matrix(0, nrow = 78*265, ncol = 1)

nday = 265
temp = 1 
for(i in 1:nstation){
  print(i)
  yhat = apply(predict_grid_covdynglg_matern_2[(0.5*M):M,i,],2, median)
  predicted_covdynglg_matern[temp:(i*nday),1] = yhat
  temp = i*nday + 1
}

temp = 1 
for(i in 1:nstation){
  print(i)
  yhat = apply(predict_grid_covdynglg[(0.5*M):M,i,],2, median)
  predicted_covdynglg[temp:(i*nday),1] = yhat
  temp = i*nday + 1
}

temp = 1 
M=10e3
for(i in 1:nday){
  print(i)
  yhat = apply(gaussian_grid$Z.loc.post[(0.5*M):M,,i],2, median)
  predicted_gaussian[temp:(i*nstation),1] = yhat
  temp = i*nstation + 1
}

ggplot() + geom_sf(df) +
  geom_point(data = info.sites4, aes(x = long, y = lat), col = 2,
             shape = 2) 


# Out-of-sample locations for spatial predictions
time_out = 266:275

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


oos_data %>% arrange(month, day)%>% 
  slice_head(n=nstation*nday) %>% arrange(code)

## Deviations between predictions of the model

all = cbind(oos_data %>% group_by(code) %>% slice_head(n=nday), 
            diff =  abs(predicted_covdynglg_matern-predicted_covdynglg))

all = cbind(oos_data %>% group_by(code) %>% slice_head(n=nday), 
            diff =  abs(predicted_covdynglg_matern-predicted_covdynglg))


all <- all %>%arrange(long)%>% group_by(code) %>% summarise(ave_diff = mean(diff), 
                                                            var_diff = var(diff),
                                                            long = long, lat= lat)


df <- map_data("world") |> 
  filter(region == "UK") #|>
  #filter(subregion %in% c("Wales", "Great Britain"))

uk <- st_as_sf(df, coords = c("long", "lat"), crs = 4326) |>
  group_by(group) |> 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

d_sf_points <- st_as_sf(all, coords = c("long", "lat"), crs = 4326)

vor <- st_voronoi(st_union(d_sf_points))
vor_poly <- st_collection_extract(vor) |>
  st_as_sf()

grid_intersect <- st_join(vor_poly, d_sf_points)

poly_intersect <- st_intersection(grid_intersect, uk)

theme_set(theme_bw(base_size = 16))

ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = ave_diff)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat), 
             shape = 2) + labs(fill = "Average of \nabsolute \ndeviations")+ 
  scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)))
 # geom_sf(data = d_sf_points, size = 1)

ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = var_diff)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat),
             shape = 2) + labs(fill = "Variance of \nabsolute \ndeviations")+ 
  scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)))

## Predictions at specific time-points

# d = 30,180, 260
d=180
all = cbind(long = unique(oos_data$long), lat = unique(oos_data$lat), 
            exp = apply(predict_grid_covdynglg[(0.5*M):M,,d],2,median),
            matern = apply(predict_grid_covdynglg_matern_2[(0.5*M):M,,d],2,median))


all <- data.frame(all) %>%arrange(long)%>%  summarise(exp = exp, matern = matern,
                                                            #ave_covdynglg = mean(covdynglg),
                                                            #ave_matern = mean(matern),
                                                            long = long, lat= lat)


df <- map_data("world") |> 
  filter(region == "UK") #|>
#filter(subregion %in% c("Wales", "Great Britain"))

uk <- st_as_sf(df, coords = c("long", "lat"), crs = 4326) |>
  group_by(group) |> 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

d_sf_points <- st_as_sf(all, coords = c("long", "lat"), crs = 4326)

vor <- st_voronoi(st_union(d_sf_points))
vor_poly <- st_collection_extract(vor) |>
  st_as_sf()

grid_intersect <- st_join(vor_poly, d_sf_points)

poly_intersect <- st_intersection(grid_intersect, uk)

theme_set(theme_bw(base_size = 16))

ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = exp)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat), col = 1, 
             shape = 2) +labs(fill = "Max ozone", caption = "CovDynGLG-Exp") + 
scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)),limits = c(35,110))#+ scale_fill_continuous(limits=c(41,93))


ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = matern)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat), col = 1,
             shape = 2) + labs(fill = "Max ozone", caption = "CovDynGLG-Matern")+
scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)),limits = c(35,110))#+ scale_fill_continuous(limits=c(41,93))



### Average of predictions over time

all = cbind(oos_data %>% group_by(code) %>% slice_head(n=nday),
            covdynglg = predicted_covdynglg, matern = predicted_covdynglg_matern)


all <- all %>%arrange(long)%>% group_by(code) %>% summarise(ave_covdynglg = mean(covdynglg),
                                                            ave_matern = mean(matern),
                                                            long = long, lat= lat)


df <- map_data("world") |> 
  filter(region == "UK") #|>
#filter(subregion %in% c("Wales", "Great Britain"))

uk <- st_as_sf(df, coords = c("long", "lat"), crs = 4326) |>
  group_by(group) |> 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

d_sf_points <- st_as_sf(all, coords = c("long", "lat"), crs = 4326)

vor <- st_voronoi(st_union(d_sf_points))
vor_poly <- st_collection_extract(vor) |>
  st_as_sf()

grid_intersect <- st_join(vor_poly, d_sf_points)

poly_intersect <- st_intersection(grid_intersect, uk)

theme_set(theme_bw(base_size = 16))


ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = ave_covdynglg)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat), 
             shape = 2) + labs(fill = "Max ozone", caption = "CovDynGLG-Exp")+
  scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)),limits = c(41,94))#+ scale_fill_continuous(limits=c(41,93))
# geom_sf(data = d_sf_points, size = 1)

ggplot() +
  geom_sf(poly_intersect, colour = NA, mapping = aes(fill = ave_matern)) +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(data = info.sites4, aes(x = long, y = lat),
             shape = 2) + labs(fill = "Max ozone", caption = "CovDynGLG-Matern")+
  scale_fill_gradientn(colours = rev(colorspace::heat_hcl(7)),limits = c(41,94))#+ scale_fill_continuous(limits=c(41,93))


## Model diagnostics

median(covdynglg_matern_grid$theta[(10e3:20e3),10])
quantile(covdynglg_matern_grid$theta[(10e3:20e3),10], prob = 0.025)
quantile(covdynglg_matern_grid$theta[(10e3:20e3),10], prob = 0.975)


### Compare coefficients of the model

M1=20e3
days = seq(as.Date("2017-3-1"), as.Date("2017-11-30"), by = "days")[1:265]
indx=1:265
data.frame(rbind(cbind(est = apply(covdynglg_grid$d.post$d0[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_grid$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_grid$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
           cbind(est = apply(covdynglg_matern_grid_2$d.post$d0[(0.5*M1):M1, ],2,median),
                 lower = apply(covdynglg_matern_grid_2$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .025),
                 upper = apply(covdynglg_matern_grid_2$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .975),
                 model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper, fill =model),alpha=.2)+
  labs(y=bquote(theta["0t"]), caption = "Intercept", fill = "Model", col = "Model") +
  theme(legend.position = c(0.85,.85))




data.frame(rbind(cbind(est = apply(covdynglg_grid$d.post$d1[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_grid$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_grid$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_grid_2$d.post$d1[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_matern_grid_2$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_grid_2$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper, fill =model),alpha=.2)+
  labs(y=bquote(theta["1t"]), caption = "Latitude") +
  theme(legend.position = "none")





data.frame(rbind(cbind(est = apply(covdynglg_grid$d.post$d2[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_grid$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_grid$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_grid_2$d.post$d2[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_matern_grid_2$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_grid_2$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper, fill =model),alpha=.2)+
  labs(y=bquote(theta["2t"]), caption = "Longitude") +
  theme(legend.position = "none")




data.frame(rbind(cbind(est = apply(covdynglg_grid$d.post$d3[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_grid$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_grid$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_grid_2$d.post$d3[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_matern_grid_2$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_grid_2$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper,fill = model),alpha=.2)+
  labs(y=bquote(theta["3t"]), caption = "Max Temperature") +
  theme(legend.position = "none")


## Plot estimated spatial precision  mixing variable
data.frame(rbind(cbind(lambda1 = apply(covdynglg_matern_grid_2$lambda1.post[(0.5*M1):M1,],2,median), model  =  "Matern", info.sites4),
           cbind(lambda1 = apply(covdynglg_grid$lambda1.post[(0.5*M1):M1,],2,median), model ="Exp", info.sites4))) %>%
  ggplot() +
  geom_sf(fill = "transparent", data = uk) +
  geom_point(aes(x=long,y=lat, col = lambda1, size = lambda1))+
  facet_wrap(~model)+
  theme(legend.position = "bottom") +
  scale_color_gradientn(colours = colorspace::heat_hcl(7))


