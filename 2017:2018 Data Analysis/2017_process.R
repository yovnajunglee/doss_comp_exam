## Process results for 2017 data  
library(ggplot2)
library(dplyr)
## Estimated coefficients to compare models

theme_set(theme_bw(base_size = 16))

days = seq(as.Date("2017-3-1"), as.Date("2017-11-30"), by = "days")[1:265]

# Need to show the estimated coefficients for covdynglg

## Estimated coefficients to compare models

M1=20e3
data.frame(rbind(cbind(est = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg5$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est)) +
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+
  labs(y=bquote(theta["0t"]), caption = "Intercept") 


data.frame(rbind(cbind(est = apply(covdynglg$d.post$d1[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg5$d.post$d1[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est))  + geom_line(aes(y = est)) +
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+
  labs(y=bquote(theta["1t"]), caption = "Latitude") 


data.frame(rbind(cbind(est = apply(covdynglg$d.post$d2[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg$d.post$d2[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est))  + geom_line(aes(y = est)) + 
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+
  labs(y=bquote(theta["2t"]), caption = "Longitude") 


data.frame(rbind(cbind(est = apply(covdynglg$d.post$d3[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg$d.post$d3[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est))  + geom_line(aes(y = est)) +
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+
  labs(y=bquote(theta["2t"]), caption = "Temperature") 




data.frame(rbind(cbind(est = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est))  + geom_line(aes(y = est)) + 
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+
  labs(y=bquote(theta["3t"]), caption = "Wind speed") 


## Temporal effects in precision model
data.frame(rbind(cbind(est = apply(covdynglg$dd.post[(0.5*M1):M1, 1,],2,median),
                       lower = apply(covdynglg$dd.post[(0.5*M1):M1, 1,],2,quantile, prob = .025),
                       upper = apply(covdynglg$dd.post[(0.5*M1):M1, 1,],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est)) + 
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+ geom_hline(yintercept = 0, linetype = 2)+
  labs(y="Temperature")

## Temporal effects in precision model
data.frame(rbind(cbind(est = apply(covdynglg$dd.post[(0.5*M1):M1, 2,],2,median),
                       lower = apply(covdynglg$dd.post[(0.5*M1):M1, 2,],2,quantile, prob = .025),
                       upper = apply(covdynglg$dd.post[(0.5*M1):M1, 2,],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = as.Date(days)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est)) + 
  geom_ribbon(aes(ymin=lower, ymax = upper),alpha=.2)+geom_hline(yintercept = 0, linetype = 2)+
  labs(y="Wind speed") + coord_cartesian(ylim = c(-.5, 1))



## Show forecasts for stations for Covdynglg and Gaussian
## at stations and months as done in paper
#indx = 1:31
#indx = 94:122
indx = 154:184
M1=1e4
M2=20e3

indout=c(6,38)

data.frame(rbind(cbind(lower = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,quantile, prob = .025),
                       fit = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,median),
                       upper = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,quantile, prob = .975),
                       model = "Gaussian", Observed = vecz[indx ,indout[1]], x = indx ),
                 cbind(lower = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2,quantile, prob = .025),
                       fit = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2, median),
                       upper = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG", Observed = vecz[indx ,indout[1]], x = indx ))) %>%
  mutate(lower = as.numeric(lower), upper = as.numeric(upper), Observed = as.numeric(Observed), fit = as.numeric(fit),
         x = as.numeric(x), Day =  rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = Observed )) + geom_point()+
  geom_line(aes(y = lower, col = model)) +
  geom_line(aes(y=upper, col = model)) +
   labs(y="Max ozone", col = "Model") + 
  #theme(legend.position = c(0.85,.15))
theme(legend.position = "none")

data.frame(rbind(cbind(lower = apply(gaussian$Z.loc.post[(0.5*M1):M1,2,indx ],2,quantile, prob = .025),
                       fit = apply(gaussian$Z.loc.post[(0.5*M1):M1,2,indx ],2,median),
                       upper = apply(gaussian$Z.loc.post[(0.5*M1):M1,2,indx ],2,quantile, prob = .975),
                       model = "Gaussian", Observed = vecz[indx ,indout[2]], x = indx ),
                 cbind(lower = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2,quantile, prob = .025),
                       fit = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2, median),
                       upper = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG", Observed = vecz[indx ,indout[2]], x = indx ))) %>%
  mutate(lower = as.numeric(lower), upper = as.numeric(upper), Observed = as.numeric(Observed), fit = as.numeric(fit),
         x = as.numeric(x), Day =  rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = Observed )) + geom_point()+
  geom_line(aes(y = lower, col = model)) +
  geom_line(aes(y=upper, col = model)) +
  labs(y="Max ozone", col = "Model") + theme(legend.position = "none")


## Posterior predictive sd
data.frame(rbind(cbind(lower = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,quantile, prob = .025),
                       fit = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,sd),
                       upper = apply(gaussian$Z.loc.post[(0.5*M1):M1,1,indx ],2,quantile, prob = .975),
                       model = "Gaussian", Observed = vecz[indx ,indout[1]], x = indx ),
                 cbind(lower = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2,quantile, prob = .025),
                       fit = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2, sd),
                       upper = apply(covdynglg$Z.loc.post[(0.5*M2):M2,1,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG", Observed = vecz[indx ,indout[1]], x = indx ))) %>%
  mutate(lower = as.numeric(lower), upper = as.numeric(upper), Observed = as.numeric(Observed), fit = as.numeric(fit),
         x = as.numeric(x), Day =  rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = Observed )) + #geom_point()+
  geom_line(aes(y = fit, col = model)) +
  #geom_line(aes(y=upper, col = model)) +
  labs(y="Posterior predictive sd", col = "Model") + 
  theme(legend.position = c(0.85,.85))
  #theme(legend.position = "none")

M1=M2=20e3
indx=1:265
data.frame(rbind(cbind(lower = apply(covdynglg_matern_nu_est$Z.loc.post[(0.5*M1):M1,2,indx ],2,quantile, prob = .025),
                       fit = apply(covdynglg_matern_nu_est$Z.loc.post[(0.5*M1):M1,2,indx ],2,sd),
                       upper = apply(covdynglg_matern_nu_est$Z.loc.post[(0.5*M1):M1,2,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG-Matern", Observed = vecz[indx ,indout[2]], x = indx ),
                 cbind(lower = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2,quantile, prob = .025),
                       fit = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2, sd),
                       upper = apply(covdynglg$Z.loc.post[(0.5*M2):M2,2,indx ],2,quantile, prob = .975),
                       model = "CovDynGLG", Observed = vecz[indx ,indout[2]], x = indx ))) %>%
  mutate(lower = as.numeric(lower), upper = as.numeric(upper), Observed = as.numeric(Observed), fit = as.numeric(fit),
         x = as.numeric(x), Day =  rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = Observed )) + #geom_point()+
  geom_line(aes(y = fit, col = model)) +
  #geom_line(aes(y=upper, col = model)) +
  labs(y="Posterior predictive sd", col = "Model") + 
  theme(legend.position = c(0.75,.85))


## Estimates of covariates in spatial precision model for full model

apply(full$beta.post[(0.5*M1):M1,],2,median)
apply(full$beta.post[(0.5*M1):M1,],2,quantile, prob = .025)
apply(full$beta.post[(0.5*M1):M1,],2,quantile, prob = .975)


### Compare Matern vs Exponential models

data.frame(rbind(cbind(est = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d0[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_matern_nu_est$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_nu_est$d.post$d0[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper, fill =model),alpha=.2)+
  labs(y=bquote(theta["0t"]), caption = "Intercept") 



data.frame(rbind(cbind(est = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG", x = 1:265),
                 cbind(est = apply(covdynglg_matern_nu_est$d.post$d4[(0.5*M1):M1, ],2,median),
                       lower = apply(covdynglg_matern_nu_est$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .025),
                       upper = apply(covdynglg_matern_nu_est$d.post$d4[(0.5*M1):M1, ],2,quantile, prob = .975),
                       model = "CovDynGLG Matern", x  = 1:265)))%>%
  mutate(est = as.numeric(est), lower = as.numeric(lower), upper = as.numeric(upper),
         Day = rep(as.Date(days)[indx],2)) %>%
  ggplot(aes(x = Day, y = est)) + geom_line(aes(y = est, col = model)) +
  geom_ribbon(aes(ymin=lower, ymax = upper, fill =model),alpha=.2)+
  labs(y=bquote(theta["0t"]), caption = "Intercept") 


median(covdynglg_matern_nu_est$theta[10e3:20e3,9])
median(covdynglg_matern_nu_est$theta[10e3:20e3,10])
quantile(covdynglg_matern_nu_est$theta[10e3:20e3,10], prob = .025)
quantile(covdynglg_matern_nu_est$theta[10e3:20e3,10], prob = .975)

## Estimated lambda 2
data.frame(rbind(cbind(lambda2 = apply(covdynglg$lambda2.post[(0.5*M2):M2,],2, median), model  = "CovDynGLG"),
           cbind(lambda2 = apply(covdynglg_matern_nu_est$lambda2.post[(0.5*M2):M2,],2, median), model  = "CovDynGLG-Matern")))%>%
  mutate(Day = rep(as.Date(days),2), lambda2 = as.numeric(lambda2))%>%
  ggplot(aes(y = lambda2)) + geom_line(aes(x=Day,colour = model))+
  geom_hline(yintercept = 1, linetype =2, col = "red") +
  labs(y = bquote(lambda["2t"]), x = "Day")
