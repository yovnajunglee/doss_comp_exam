
## Fit Full model
## Code adapted from Fonseca et al. (2023) 
## https://github.com/thaiscofonseca/DynGLG


run_full_model <- function(M=10){
  
  theta= matrix(0,M,9)  ## sig2, phi, alpha, nu, like, like2
  # MCMC samples for coefficients in level 1
  d0.post<- matrix(NA, M, ncol=J)
  d1.post<- matrix(NA, M, ncol=J)
  d2.post<- matrix(NA, M, ncol=J)
  d3.post<- matrix(NA, M, ncol=J)
  d4.post<- matrix(NA, M, ncol=J)
  # MCMC samples for mixing process
  # J = no. of days
  # Inew = no. of stations
  lambda.post<- array(data=NA, c(M,Inew,J))
  mu.post<- matrix(NA, M, ncol=J) # ???
  # MCMC samples for spatial mixing process
  lambda1.post<- array(data=NA, c(M,Inew))
  # MCMC samples for temporal mixing process
  lambda2.post<- array(data=NA, c(M,J))
  ##dd.post<- array(data=NA, c(M,1,J)) #### MUDAR AQUI
  
  ## Variance tuning for RW - MH
  ra = .05;  rnu1= .7; rnu2 = 0.5 ; ral=0.04 ; raa = .2; 
  rs2=0.1 ; rt2=0.15
  
  rlamb = c(0.02,0.02,0.02, 0.02)*.3
  rlamb2 = rep(0.09,10) 
  
  cont1 = 0; cont2=0;cont3=0; cont4=rep(0,4); cont3b = 0;  ## por alpha and phi
  cont4b<- rep(0,10)
  
  
  
  #------------------------------------------------------------------
  ### initial values
  sig2k = mean(apply(z,1,var));  # Sigma2 in eq.(1) # Gaussian noise 
  #  0.001;   # Nugget effect
  tau2k = 0.1*sig2k
  nuk=0.15; ak = 0.5;alphak= 1
  C1<- CauchyCov(distS, ak, alphak)
  
  
  nu1k = 0.1 # eq. (3)
  #nu1k = 0
  loglamb1k = rep(0,Inew)
  nu2k = 0.1 # eq. (4a)
  #nu2k = 0
  mutk = rep(0,J) ## ??
  loglamb2k = rep(0,J)
  
  # Covariates in the spatial mixing process
  X.lamb1 =  cbind(lat,long)
  X.lamb1.out = cbind(lat.out,long.out)
  #X.lamb1 <- cbind((rep(1,Inew)))
  betak = matrix(0,1,dim(X.lamb1)[2]) 
  beta.post = matrix(NA,M,dim(X.lamb1)[2])
  
  
  loglambk = matrix(0,J,Inew)
  # Log likelihood for log(lambda) = log(lambda1) + log(lambda_2)
  for (i in 1:Inew){ for (j in 1:J){ loglambk[j,i] <- loglamb1k[i]+loglamb2k[j] }}
  
  
  
  #-----------------------------------------------------------------
  #### FFBS information
  #### for theta.t
  Ftt = NULL
  #List where each item is covariates for each timepoint
  for (j in 1:J){
    Ftt[[j]] = t(cbind(rep(1,Inew), Plat, Plong, Ptemp[j,],Pws[j,])) ##  n x r
  }
  R = dim(Ftt[[1]])[1]
  Gt<- diag(R)  ### p+1 covariates    ### n x n 
  C0<- 100*diag(R)
  m0<- rep(0,R)
  discount<- c(.99)
  ### for mu.t
  #m1<- 0 ; C11<- 100 ; Ft1<- t(rep(1,Inew)); Gt1<- 1  
  d2<- c(.99)  ## n x n
  
  # Design matrix containing mean temperature and wind
  # for each day over the sites
  Ft1 = NULL
  Ft1 = t( cbind(apply(Ptemp,1,mean),apply(Pws,1,mean))) ##  n x r

  Ft1.pred = t(cbind(apply(Ptemp.p,1,mean),apply(Pws.p,1,mean)))
  #Ft1 = t(cbind(rep(1,J))) ##  n x 1 (intercept only)
  
  # I believe is eta_t from eq. (4a) Yes = eta
  dd.post<- array(data=NA, c(M,dim(Ft1)[1],J))
  
  R1 = dim(Ft1)[1]
  Gt1<- diag(R1)  ### p+1 covariates    ### n x n 
  C11<- 100*diag(R1)
  m1<- rep(0,R1)
  ##################################################
  In<-diag(1, Inew)  ## matriz diagonal de 1's
  
  # posterior samples of theta
  h=10
  theta.h.post =  array(data=NA, c(M,nrow(Ftt[[1]]),h))
  eta.h.post =  array(data=NA, c(M,dim(Ft1)[1],h))
  Z.pred.post = array(data=NA, c(M,Inew,h))
  Z.loc.post = array(data=NA, c(M,length(Plat.out),J))
  #*********************************** MCMC *****************************/
  #********************* Non-Gaussian *****************************/
  ######### criandos regioes (blocos) para lambda
  
  ind = 1:Inew
  region3 = ind[which((coords[,1]< (52.5))&(coords[,2]< (-1.8)))] 
  region2 = ind[which((coords[,1]>= (52.5))&(coords[,2]>= (-1.8)))] 
  region1 = ind[which((coords[,1]>= (52.5))&(coords[,2]< (-1.8)))] 
  region4 = ind[which((coords[,1]< (52.5))&(coords[,2]>= (-1.8)))]
  
  indr = c(region1,region2,region3,region4)
  n<-NULL
  n[1] = length(region1) 
  n[2] = length(region2) 
  n[3] = length(region3) 
  n[4] = length(region4) 
  nr = c(0,cumsum(n))
  
  indrt = 1:J
  nt = c(rep(30,8),J-sum(rep(30,8)))
  nrt = c(0,cumsum(nt))

  epsk = z*0
  d0k=d1k=d2k=d3k=d4k=rep(0,J)
  
  # obtain log likelihood for eq. (1)
  lverok = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k)   ### d0k,d1k,d2k,d3k
  
  # Log-likelihood for log(lambda_1)
  # Regression for spatial mixing process
  lverok1 = LogVeroLogLamb1reg(loglamb1k,distS,ak,alphak,nu1k,betak,X.lamb1)
  # loglikelihood for lambda_2 temporal 
  # mutk = F2*eta
  lverok2 = LogVeroLogLamb2(loglamb2k,nu2k,mutk)
  
  a1k = 10
  alpha1k = 0.5
  sig2k = var(as.vector(z))
  
  
  
  system.time(for (k in 1:M){
    
    Cc = CauchyCov(distS,ak,alphak)   
    Cinv = solve(Cc)
    ## FFBS for deltas (each time t)
    ### TROCAR AQUI O FATOR DE DESCONTO
    ## Run FBBS to obtain theta_t (note conditionally Gaussian)
    ffbs.theta <- FFBS(z-epsk,Cc,discount[1],m0,C0,Ftt,Gt, loglamb1k,loglamb2k,sig2k,tau2k, In) 
    dt <- ffbs.theta$d
    d0.post[k,]<- dt[1,1,1:J] ;   d1.post[k,]<- dt[2,1,1:J] ;   
    d2.post[k,]<- dt[3,1,1:J] ;   d3.post[k,]<- dt[4,1,1:J];
    d4.post[k,]<- dt[5,1,1:J]
    
    Sinv = diag(as.vector(sqrt(exp(loglamb1k))))%*%Cinv%*%diag(as.vector(sqrt(exp(loglamb1k))))
    for (j in 1:J){
      mu = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +  dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      aux = solve(exp(loglamb2k[j])*Sinv/sig2k+In/tau2k)
      epsk[j,] = mvrnorm(1,aux%*%(exp(loglamb2k[j])*Sinv/sig2k)%*%(z[j,]-mu),aux) #what is eps?
    }
    # store theta_t+h
    #  print(dim( theta.h.post[k,,]));  print(dim( ffbs.theta$theta.h))
    theta.h.post[k,,] <- ffbs.theta$theta.h 
    
    
    ### FFBS for mus (each time t)  - lambda_t
    
    ## Now mu = F_2*eta
    V1<- nu2k*1
    z.lk<- matrix(loglamb2k + (nu2k/2),J,1)
    
    # eta
    ffbs.mu = FFBS.mu(z.lk, V1, d2[1], m1, C11, Ft1, Gt1)
    aux = ffbs.mu$d
    # F_2*eta
    mu.post[k,] <- apply(aux[1:R1,1,1:J]*Ft1,2,sum)
    
    # eta
    dd.post[k,,] = aux[,1,]
    #print(dim(eta.h.post[k,,]))
    #print(dim(ffbs.mu$eta.h ))
    eta.h.post[k,,] <- ffbs.mu$eta.h 
    
    # Update current value
    d0k<- d0.post[k,] ; d1k<- d1.post[k,] ; d2k<- d2.post[k,]; d3k<-d3.post[k,]; d4k<-d4.post[k,]
    mutk<- mu.post[k,]
    
    # MH Sampling to obtain lambda1 
    # Proposal candidate for loglambda1
    loglamb1prop = rep(NA, Inew)
    ###################### generate lambda(s,t)
    for (j in 1:4){
      
      region = indr[(nr[j]+1):nr[j+1]]
      loglamb1prop = loglamb1k;
      # Proposal value for loglamb1
      # Random walk MH
      loglamb1prop[region] = loglamb1k[region] + rlamb[j]*rnorm(length(region));
      
      # log-likelihood at current value
      lverok = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k)   ### d0k,d1k,d2k,d3k
      lverok1 = LogVeroLogLamb1reg(loglamb1k,distS,a1k,alpha1k,nu1k,betak,X.lamb1)
      # log likelihood at proposal value
      lveroprop = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1prop,loglamb2k,tau2k);    ### d0k,d1k,d2k,d3k
      lveroprop1 = LogVeroLogLamb1reg(loglamb1prop,distS,a1k,alpha1k,nu1k,betak,X.lamb1)
      
      auxprop = lveroprop+lveroprop1 ;
      auxk = lverok+lverok1 ;
      # acceptance probability
      ratio = auxprop-auxk; test = runif(1);
      
      if (ratio>log(test)) {
        loglamb1k[region] = loglamb1prop[region];
        lverok1 = lveroprop1; lverok = lveroprop;
        cont4[j]<- cont4[j]+1
      }
      
      lambda1.post[k,]<- exp(loglamb1k)   ## array
      
      
    }
    
    ##### covariates in lambda1
    ### Change this if we dont need betak 
    
    y.lamb1 = loglamb1k+nu1k/2
    S.lamb1 = nu1k*exp(-distS/a1k)
    P.lamb1 = solve(S.lamb1)
    auxk = solve(t(X.lamb1)%*%P.lamb1%*%X.lamb1)
    bet.hat = auxk%*%t(X.lamb1)%*%P.lamb1%*%y.lamb1
    ## Gibbs sampling for beta
    #betak=0*bet.hat
    betak = rmvnorm(1,mean=bet.hat,sigma=auxk)
    #betak=0
    
    # RW MH for lambda2
    loglamb2prop = rep(NA, J)
    ###################### generate lambda(s,t)
    for (j in 1:9){
      
      regiont = indrt[(nrt[j]+1):nrt[j+1]]
      loglamb2prop = loglamb2k;
      # Current value + Gaussian noise
      loglamb2prop[regiont] = loglamb2k[regiont] + rlamb2[j]*rnorm(length(regiont));
      
      lverok = LogVeroNGt(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k,regiont)   ### d0k,d1k,d2k,d3k
      lverok2 = LogVeroLogLamb2(loglamb2k[regiont],nu2k,mutk[regiont])
      lveroprop = LogVeroNGt(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2prop,tau2k,regiont);    ### d0k,d1k,d2k,d3k
      lveroprop2 = LogVeroLogLamb2(loglamb2prop[regiont],nu2k,mutk[regiont])
      
      auxprop = lveroprop+lveroprop2 ;
      auxk = lverok+lverok2 ;
      ratio = auxprop-auxk; test = runif(1);
      
      if (ratio>log(test)) {
        loglamb2k[regiont] = loglamb2prop[regiont];
        cont4b[j]<- cont4b[j]+1
      }
      
      lambda2.post[k,]<- exp(loglamb2k)   ## array
      
      
    }
    ### gibbs for sig2
    
    bk = 0
    Cinv = solve(CauchyCov(distS,ak,alphak))
    Sinv = diag(as.vector(sqrt(exp(loglamb1k))))%*%Cinv%*%diag(as.vector(sqrt(exp(loglamb1k))))
    for (j in 1:J){
      mu = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +  dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      media = z[j,]-mu-epsk[j,]
      bk = bk + t(media)%*%(exp(loglamb2k[j])*Sinv)%*%media
    }
    sig2k = 1/rgamma(1,(J*Inew)/2+.0001,bk/2 + 0.0001)
    
    ### gibbs for tau2
    
    bk = 0
    for (j in 1:J){bk = bk + t(epsk[j,])%*%epsk[j,]}
    tau2k = 1/rgamma(1,(J*Inew)/2+0.0001,bk/2 + 0.001)
    

    ### M-H for phi and alpha (do both at the same time)
    
    aprop = exp(log(ak)+ra*rnorm(1)); # Proposed value for phi = scale and shape param in Cauchy corr. function in obs. eq.
    alphaAux = log(alphak/(2-alphak))+ral*rnorm(1); 
    alphaprop = 2*exp(alphaAux)/(1+exp(alphaAux));# Proposed value for alpha
    
    
    lverok = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k)
    lveroprop = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,aprop,alphaprop,loglamb1k,loglamb2k,tau2k)
    
    # distribution at proposal - c(Phi,alpha)
    auxprop = lveroprop+log(dgamma(aprop,1, 2.3/median(distS))) + log(aprop)+ log(alphaprop*(2-alphaprop))    #
    auxk = lverok+log(dgamma(ak,1, 2.3/median(distS)))  + log(ak) +   log(alphak*(2-alphak))  ##
    
    ratio = auxprop-auxk;   test = runif(1);
    if (ratio>log(test)) {
      ak = aprop;
      alphak = alphaprop;
      lverok = lveroprop;
      cont2 = cont2 + 1;
    }
    
    ### Metropolis-Hastings for nu and gamma in exp. corr. function in spatial mixing process.
    a1prop = exp(log(a1k)+raa*rnorm(1)); ## a1 = gamma in covariance
    nu1prop = exp(log(nu1k)+rnu1*rnorm(1)); ## nu1
    lverok1 = LogVeroLogLamb1reg(loglamb1k,distS,a1k,alpha1k,nu1k,betak,X.lamb1)
    lveroprop1 = LogVeroLogLamb1reg(loglamb1k,distS,a1prop, alpha1k, nu1prop,betak,X.lamb1)
    
    auxprop = lveroprop1+log(nu1prop)+  log(dgamma(nu1prop,1,5))+log(a1prop)+  log(dgamma(a1prop,1,5))    #log(dgig(nuprop, 0, 0.75,6))  #
    auxk = lverok1+log(nu1k)+   log(dgamma(nu1k,1,5)) +log(a1k)+   log(dgamma(a1k,1,5))
    ratio = auxprop-auxk;   test = runif(1);
    if (ratio>log(test)) {
      nu1k = nu1prop;
      a1k = a1prop;
      lverok1 = lveroprop1;
      cont3 = cont3 + 1;
    }
    
    ####### nu2
    
    nu2prop = exp(log(nu2k)+rnu2*rnorm(1));
    lverok2 = LogVeroLogLamb2(loglamb2k, nu2k,mutk)
    lveroprop2 = LogVeroLogLamb2(loglamb2k, nu2prop,mutk)
    
    auxprop = lveroprop2+log(nu2prop)+  log(dgamma(nu2prop,1,5))    #log(dgig(nuprop, 0, 0.75,6))  #
    auxk = lverok2+log(nu2k)+   log(dgamma(nu2k,1,5)) #;log(dgig(nuk, 0, 0.75,6)) #
    ratio = auxprop-auxk;   test = runif(1);
    if (ratio>log(test)) {
      nu2k = nu2prop;
      lverok2 = lveroprop2;
      cont3b = cont3b + 1;
    }
    
    ## Forecasts at obs. locations
    
    Cc = CauchyCov(distS,ak,alphak)   
    Cinv = solve(Cc)
    lambda2.pred <- exp(mvrnorm(1,apply(eta.h.post[k,,]*Ft1.pred,2,sum) - .5*nu2k, nu2k*diag(h) ))
    ## Generate Z_t_pred [for avail. locations]
    Sinv = diag(as.vector(sqrt(exp(loglamb1k))))%*%Cinv%*%diag(as.vector(sqrt(exp(loglamb1k))))
    for (hh in 1:h){
      mu = theta.h.post[k,1,hh] +  theta.h.post[k,2,hh]*Plat +  theta.h.post[k,3,hh]*Plong +  
        theta.h.post[k,4,hh]*Ptemp.p[hh,] +  theta.h.post[k,5,hh]*Pws.p[hh,]
      #print(length(theta.h.post))
      
      aux = solve((lambda2.pred[hh])*Sinv/sig2k +In/tau2k)
      #print(dim(aux))
      
      Z.pred.post[k,,hh] = mvrnorm(1,mu,aux)
    }
    
    ## Spatial interpolation at unobserved sites
    
    ## a) Spatial mixing process
    
    Wlamb.out = exp(-distS.out/a1k)
    Winv.out = solve(Wlamb.out)
    Wlamb.p.out = exp(-distS.tot/a1k)[1:Inew, (Inew+1):(Inew + length(Plat.out))]
    a = loglamb1k + 0.5*nu1k - (X.lamb1%*%t(betak))
    Wlamb= exp(-distS/a1k)  ## covariance matrix for lambda
    Winv<- solve(Wlamb)
    mu.loglamb.out = -nu1k/2*rep(1, dim(distS.out)[1]) +   (X.lamb1.out%*%t(betak)) + 
      t(Wlamb.p.out)%*% t(Winv)%*%a
    Sig.loglamb.out = nu1k*( Wlamb.out - (t(Wlamb.p.out)%*%Winv%*%(Wlamb.p.out)))
    loglamb1.out.k = mvrnorm(1, mu.loglamb.out, Sig.loglamb.out)
    
    
    ## b) Z - observation process
    
    Cc = CauchyCov(distS,ak,alphak)
    Cinv = solve(Cc)
    S = diag(as.vector(1/sqrt(exp(loglamb1k))))%*%Cc%*%diag(as.vector(1/sqrt(exp(loglamb1k))))
    
    Cc.out =  CauchyCov(distS.out,ak,alphak)
    Cinv.out = solve(Cc.out)
    S.out = (diag(as.vector(1/sqrt(exp(loglamb1.out.k))))%*%Cc.out%*%
               diag(as.vector(1/sqrt(exp(loglamb1.out.k)))))
    
    
    Cc.p.out = CauchyCov(distS.tot,ak,alphak)
    S.p.out = diag(as.vector(1/sqrt(exp(c(loglamb1k, loglamb1.out.k)))))%*%
      Cc.p.out%*%diag(as.vector(1/sqrt(exp(c(loglamb1k, loglamb1.out.k)))))
    
    
    
    Z.pred.k <- matrix(NA, nrow = J, ncol = length(Plat.out))
    
    for (j in 1:J){
      mu.obs = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +
        dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      mu.pred = dt[1,1,j] +  dt[2,1,j]*Plat.out +  dt[3,1,j]*Plong.out +
        dt[4,1,j]*Ptemp.out[j,]+  dt[5,1,j]*Pws.out[j,]
      b = z[j,] - mu.obs
      
      aux.out = sig2k*(1/exp(loglamb2k[j]))*S.out + tau2k*diag(1, length(Plat.out))
      auxinv.out = solve(aux.out)
      aux = sig2k*(1/exp(loglamb2k[j]))*S + tau2k*diag(1, Inew)
      auxinv = solve(aux)
      
      aux.p.out = sig2k*(1/exp(loglamb2k[j]))*S.p.out + tau2k*diag(1, (Inew + length(Plat.out)))
      aux.p.out = aux.p.out[1:Inew, (Inew+1):(Inew + length(Plat.out))]
      cov.z = aux.out - (t(aux.p.out)%*%auxinv%*%(aux.p.out))
      
      mu.z.pred = mu.pred +(t(aux.p.out)%*%auxinv%*%b)
      
      Z.loc.post[k,,j] = mvrnorm(1, mu.z.pred, cov.z)
    }
    

    
    
    
    beta.post[k,] = betak;
    theta[k,1]= sig2k;
    theta[k,2] = ak;
    theta[k,3] = alphak;
    theta[k,4] = nu1k;
    theta[k,5] = tau2k;
    theta[k,6] = lverok;
    theta[k,7] = lverok2;
    theta[k,8] = nu2k;
    theta[k,9] = a1k;
    
    print(k)
    
    
  }) 
  
  return(list(theta=theta, 
              beta.post = beta.post,
              dd.post = dd.post,
              d.post = list(d0 = d0.post,
                            d1 = d1.post,
                            d2 = d2.post,
                            d3 = d3.post,
                            d4 = d4.post),
              theta.h.post = theta.h.post,
              eta.h.post = eta.h.post,
              Z.loc.post =Z.loc.post,
              Z.pred.post = Z.pred.post,
              lambda1.post = lambda1.post,
              lambda2.post = lambda2.post,
              lambda.post = lambda.post,
              mu.post = mu.post,
              cont = list(cont1, cont2, cont3, cont3b, cont4, cont4b))
         
  )
}