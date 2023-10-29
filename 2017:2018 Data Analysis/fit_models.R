### Process data
rm(list=ls())

# Load data
load("~/Data/2017data.RData") 
source("~/Desktop/Research Comps/Code/Author/GaussianModel.R")

vecz1 = vecz
indout = c(6,38) ## Test set
time_out = 266:275 # Forecast horizon
vecz = vecz1[-time_out,]
J= nrow(vecz) ;I = ncol(vecz) ### 303 por 61
z = vecz[,-indout]     
Inew = I-length(indout)  ##  locations (estimation) 

lat<- info.sites4$lat[-indout] ; long<- info.sites4$long[-indout] ; 
x1<-info.sites4$type[-indout]
Plat<- (lat - median(lat))/diff(quantile(lat, c(0.25,0.75))) 
Plong<- (long - median(long))/diff(quantile(long, c(0.25,0.75))) 
Ptemp = (t(dt.temp4[-indout,-time_out])-median(t(dt.temp4[-indout,-time_out])))/diff(quantile(dt.temp4[-indout,-time_out], c(0.25,0.75))) 
Pws = (t(dt.ws4[-indout,-time_out])-median(t(dt.ws4[-indout,-time_out])))/diff(quantile(dt.ws4[-indout,-time_out], c(0.25,0.75))) 


# Out-of-sample time-points for forecasts at observed locations 
Ptemp.p = (t(dt.temp4[-indout,time_out])-median(t(dt.temp4[-indout,time_out])))/diff(quantile(dt.temp4[-indout,time_out], c(0.25,0.75))) 
Pws.p = (t(dt.ws4[-indout,time_out])-median(t(dt.ws4[-indout,time_out])))/diff(quantile(dt.ws4[-indout,time_out], c(0.25,0.75))) 


# Out-of-sample locations for spatial predictions
#indout = test.set1
lat.out<- info.sites4$lat[indout] ; long.out<- info.sites4$long[indout] ; 
x1.out = info.sites4$type[indout]
Plat.out<- (lat.out - median(lat.out))/diff(quantile(lat.out, c(0.25,0.75))) 
Plong.out<- (long.out - median(long.out))/diff(quantile(long.out, c(0.25,0.75))) 
Ptemp.out = (t(dt.temp4[indout,-time_out])-median(t(dt.temp4[indout,-time_out])))/diff(quantile(dt.temp4[indout,-time_out], c(0.25,0.75))) 
Pws.out = (t(dt.ws4[indout,-time_out])-median(t(dt.ws4[indout,-time_out])))/diff(quantile(dt.ws4[indout,-time_out], c(0.25,0.75))) 

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
source("~/Desktop/Research Comps/Code/Author/Personal/run_gaussian.R")
gaussian <- run_gaussian_model(M=10000)

source("~/Desktop/Research Comps/Code/Author/Personal/run_covdynglg.R")
covdynglg <- run_covdynglg(M=20e3)

source("~/Desktop/Research Comps/Code/Author/Personal/run_covdynglg_matern_fields.R")
covdynglg_matern_nu_est <- run_covdynglg_matern(M=20000, nu =1.5)

source("~/Desktop/Research Comps/Code/Author/Personal/run_full.R")
full <- run_full_model(M=20000)

source("~/Desktop/Research Comps/Code/Author/Personal/run_studentt.R")
studentt = run_studentt(M=10e3)



## Model Diagnostics

M=20000
plot(covdynglg_matern_nu_est$theta[(0.5*M):M,10],type='l')
plot(covdynglg_matern_nu_est$theta[(0.5*M):M,9], type ='l')
plot(covdynglg_matern$theta[(0.5*M):M,4], type ='l')
median(covdynglg$theta[(0.5*20000):20000,9])
median(covdynglg_matern_05$theta[(.5*M):M,9])

## Plot spatial predictions

indx=1:265
data.frame(rbind(cbind(lower = apply(covdynglg$Z.loc.post[(0.5*M):M,1,indx ],2,quantile, prob = .025),
                       upper = apply(covdynglg$Z.loc.post[(0.5*M):M,1,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG", Observed = vecz[indx ,indout[1]], x = indx ),
                 cbind(lower = apply(covdynglg_matern_nu_est$Z.loc.post[(0.5*M):M,1,indx ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_nu_est$Z.loc.post[(0.5*M):M,1,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG matern", Observed = vecz[indx ,indout[1]], x = indx ),
                 cbind(lower = apply(gaussian$Z.loc.post[(0.5*10000):10000,1,indx ],2,quantile, prob = .025),
                       upper = apply(gaussian$Z.loc.post[(0.5*10000):10000,1,indx ],2,quantile, prob = .975),
                       model = "Gaussian", Observed = vecz[indx ,indout[1]], x = indx )
)) %>%
  mutate(lower = as.numeric(lower), upper = as.numeric(upper), Observed = as.numeric(Observed),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = Observed )) + geom_point()+ geom_line(aes(y = lower, col = model)) +
  geom_line(aes(y=upper, col = model)) +
  theme_bw() + labs(y="Max ozone")


## Look at estimated coefficients to compare coefficients
# 
data.frame(rbind(cbind(est = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,median),
                       model = "CovDynGLG", x = 1:265 ),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d0[(0.5*M):M, ],2,median),
                       model = "CovDynGLG Matern ", x = 1:265 )))%>%
  mutate(est = as.numeric(est),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = est)) + geom_line(aes(y = est, col = model)) + theme_bw() + labs(y="Intercept")


data.frame(rbind(cbind(est = apply(covdynglg$d.post$d1[(0.5*M1):M1, ],2,median),
                       model = "CovDynGLG", x = 1:265 ),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d1[(0.5*M):M, ],2,median),
                       model = "CovDynGLG Matern ", x = 1:265 )))%>%
  mutate(est = as.numeric(est),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = est)) + geom_line(aes(y = est, col = model)) + theme_bw() + labs(y="Latitude")

# 
data.frame(rbind(cbind(est = apply(covdynglg$d.post$d2[(0.5*M1):M1, ],2,median),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d2[(0.5*M):M, ],2,median),
                       model = "CovDynGLG Matern ", x = 1:265 )))%>%
  mutate(est = as.numeric(est),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = est)) + geom_line(aes(y = est, col = model)) + theme_bw() + labs(y="Longitude")

# 
data.frame(rbind(cbind(est = apply(covdynglg$d.post$d3[(0.5*M):M, ],2,median),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d3[(0.5*M):M, ],2,median),
                       model = "CovDynGLG Matern ", x = 1:265 )))%>%
  mutate(est = as.numeric(est),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = est)) + geom_line(aes(y = est, col = model)) + theme_bw() + labs(y="Temperature")

# 
data.frame(rbind(cbind(est = apply(covdynglg$d.post$d4[(0.5*M):M, ],2,median),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d4[(0.5*M):M, ],2,median),
                       model = "CovDynGLG Matern ", x = 1:265)))%>%
  mutate(est = as.numeric(est),
         x = as.numeric(x)) %>%
  ggplot(aes(x = x, y = est)) + geom_line(aes(y = est, col = model)) + theme_bw() + labs(y="Windspeed")

