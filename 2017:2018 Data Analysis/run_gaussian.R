
## Fit Gaussian model with nugget effect
## Code adapted from Fonseca et al. (2023) 
## https://github.com/thaiscofonseca/DynGLG

run_gaussian_model <- function(M=10000){
  
  theta= matrix(0,M,9)  
  d0.post<- matrix(NA, M, ncol=J)
  d1.post<- matrix(NA, M, ncol=J)
  d2.post<- matrix(NA, M, ncol=J)
  d3.post<- matrix(NA, M, ncol=J)
  d4.post<- matrix(NA, M, ncol=J)
  lambda.post<- array(data=NA, c(M,Inew,J))
  mu.post<- matrix(NA, M, ncol=J)
  lambda1.post<- array(data=NA, c(M,Inew))
  lambda2.post<- array(data=NA, c(M,J))
  dd.post<- array(data=NA, c(M,3,J)) 
  h=10
  
  ra = .07;  rnu1= 1.0; rnu2 = 0.5 ; ral=0.03 ; raa = .1; 
  rs2=0.1 ; rt2=0.08
  
  cont1 = 0; cont2=0;cont3=0; cont4=rep(0,4); cont3b = 0;  
  cont4b<- rep(0,10)
  
  
  
  #------------------------------------------------------------------
  ### initial values
  sig2k = mean(apply(z,1,var));  
  #tau2k = 0.1*sig2k
  tau2k = 0.001
  nuk=0.15; ak = 0.5;alphak= 1
  C1<- CauchyCov(distS, ak, alphak)
  
  
  nu1k = 0.1
  loglamb1k = rep(0,Inew)
  nu2k = 0.1
  mutk = rep(0,J)
  loglamb2k = rep(0,J)
  
  loglambk = matrix(NA,J,Inew)
  for (i in 1:Inew){ for (j in 1:J){ loglambk[j,i] <- loglamb1k[i]+loglamb2k[j] }}
  
  
  
  #-----------------------------------------------------------------
  #### FFBS information
  #### for theta.t
  Ftt = NULL
  for (j in 1:J){
    Ftt[[j]] = t(cbind(rep(1,Inew), Plat, Plong, Ptemp[j,],Pws[j,])) ##  n x r
  }
  R = dim(Ftt[[1]])[1]
  Gt<- diag(R)  ### p+1 covariates    ### n x n 
  C0<- 10*diag(R)
  m0<- rep(0,R)
  discount<- .95
  
  In<-diag(1, Inew)  ## matriz diagonal de 1's
  
  
  #*********************************** MCMC *****************************/
  
  epsk = z*0
  
  d0k=d1k=d2k=d3k=d4k=rep(0,J) #
  lverok = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k)   ### d0k,d1k,d2k,d3k
  
  a1k = 10
  alpha1k = 0.5
  sig2k = var(as.vector(z))
  

  theta.h.post =  array(data=NA, c(M,nrow(Ftt[[1]]),h))
  #eta.h.post =  array(data=NA, c(M,dim(Ft1)[1],h))
  Z.pred.post = array(data=NA, c(M,Inew,h))
  Z.loc.post = array(data=NA, c(M,length(Plat.out),J))
  
  
  system.time(for (k in 1:M){
    
    Cc = CauchyCov(distS,1,1)   
    Cinv = solve(Cc)
    
    ## FFBS for deltas (each time t)
    ffbs.theta <- FFBS.g(z-epsk,Cc,discount[1],m0, C0,Ftt,Gt, loglamb1k,loglamb2k,sig2k,tau2k, In)
    dt = ffbs.theta$d
    d0.post[k,]<- dt[1,1,1:J] ;   d1.post[k,]<- dt[2,1,1:J] ;   
    d2.post[k,]<- dt[3,1,1:J] ;   d3.post[k,]<- dt[4,1,1:J];
    d4.post[k,]<- dt[5,1,1:J]
    
    ### uncomment for general model with nugget effect  
    Sinv = Cinv
    aux = solve(Sinv/sig2k+In/tau2k)
    for (j in 1:J){
      mu = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +  dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      epsk[j,] = mvrnorm(1,aux%*%(Sinv/sig2k)%*%(z[j,]-mu),aux)
    }
    
    ### FFBS for mus (each time t)  - lambda_t
    mu.post[k,] <- rep(0,J)
    mutk<- mu.post[k,]
    
    
    ### uncomment for general model with nugget effect  
    ### gibbs for tau2
    bk = 0
    for (j in 1:J){bk = bk + t(epsk[j,])%*%epsk[j,]}
    tau2k = 1/rgamma(1,(J*Inew)/2+0.001,bk/2+0.001)
    
    ### gibbs for sig2
    
    bk = 0
    Cinv = solve(CauchyCov(distS,ak,alphak))
    for (j in 1:J){
      mu = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +  dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      media = z[j,]-mu-epsk[j,]
      bk = bk + t(media)%*%Cinv%*%media
    }
    sig2k = 1/rgamma(1,(J*Inew)/2+0.001,bk/2 + 0.001)
    
    ### M-H for phi and alpha
    
    aprop = exp(log(ak)+ra*rnorm(1));
    alphaAux = log(alphak/(2-alphak))+ral*rnorm(1);
    alphaprop = 2*exp(alphaAux)/(1+exp(alphaAux));
    
    lverok = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,ak,alphak,loglamb1k,loglamb2k,tau2k)
    lveroprop = LogVeroNG(z-epsk,distS,Plat,Plong,Ptemp,Pws,d0k,d1k,d2k,d3k,d4k,sig2k,aprop,alphaprop,loglamb1k,loglamb2k,tau2k)
    
    auxprop = lveroprop+log(dgamma(aprop,1, 2.3/median(distS))) + log(aprop)+ log(alphaprop*(2-alphaprop))    #
    auxk = lverok+log(dgamma(ak,1, 2.3/median(distS)))  + log(ak) +   log(alphak*(2-alphak))  ##
    
    ratio = auxprop-auxk;   test = runif(1);
    if (ratio>log(test)) {
      ak = aprop;
      alphak = alphaprop;
      lverok = lveroprop;
      cont2 = cont2 + 1;
    }
    
    theta.h.post[k,,] <- ffbs.theta$theta.h 
    

    
    ## Forecasts at obs. locations
    
    Cc = CauchyCov(distS,ak,alphak)   
    Cinv = solve(Cc)
    lambda2.pred <- exp(loglamb2k) 
    ## Generate Z_t_pred [for avail. locations]
    Sinv = diag(as.vector(sqrt(exp(loglamb1k))))%*%Cinv%*%diag(as.vector(sqrt(exp(loglamb1k))))
    for (hh in 1:h){
      mu = theta.h.post[k,1,hh] +  theta.h.post[k,2,hh]*Plat +  theta.h.post[k,3,hh]*Plong +  
        theta.h.post[k,4,hh]*Ptemp.p[hh,] +  theta.h.post[k,5,hh]*Pws.p[hh,]
      aux = solve((lambda2.pred[hh])*Sinv/sig2k+In/tau2k)
      Z.pred.post[k,,hh] = mvrnorm(1,mu,aux)
    }
    
    ## Spatial interpolation at unobserved sites
   
    loglamb1.out.k = rep(0, length(indout))
    
    ## b) Z - observation process
    
    Cc = CauchyCov(distS,ak,alphak)   
    Cinv = solve(Cc)
    Cc.out =  CauchyCov(distS.out,ak,alphak) 
    Cc.p.out = CauchyCov(distS.tot,ak,alphak)
    Z.pred.k <- matrix(NA, nrow = J, ncol = (length(Plat.out)))
    
    for (j in 1:J){
      mu.obs = dt[1,1,j] +  dt[2,1,j]*Plat +  dt[3,1,j]*Plong +
        dt[4,1,j]*Ptemp[j,]+  dt[5,1,j]*Pws[j,]
      mu.pred = dt[1,1,j] +  dt[2,1,j]*Plat.out +  dt[3,1,j]*Plong.out +
        dt[4,1,j]*Ptemp.out[j,]+  dt[5,1,j]*Pws.out[j,]
      b = z[j,] - mu.obs - epsk[j,]
      aux.out = sig2k*Cc.out + tau2k*diag(1, length(Plat.out))
      aux=  sig2k*Cc + tau2k*diag(1, Inew)
      aux.p.out = sig2k*Cc.p.out + tau2k*diag(1, Inew + length(Plat.out))
      aux.p.out = aux.p.out[1:Inew, (Inew+1):(Inew + length(Plat.out))] 
      
      auxinv = solve(aux)
      cov.z = aux.out - (t(aux.p.out)%*%auxinv%*%(aux.p.out))
      
      mu.z.pred = mu.pred +(t(aux.p.out)%*%auxinv%*%b)
      
      Z.loc.post[k,,j]  = mvrnorm(1, mu.z.pred, cov.z)
      
    }
    

    
    
    
    theta[k,1]= sig2k;
    theta[k,2] = ak;
    theta[k,3] = alphak;
    theta[k,4] = nu1k;
    theta[k,5] = tau2k;
    theta[k,6] = lverok;
    theta[k,7] = 0;# lverok2;
    theta[k,8] = nu2k;
    theta[k,9] = a1k;
    
    print(k)
    
    
  })
  return(list(theta=theta, 
              #beta.post = beta.post,
              dd.post = dd.post,
              d.post = list(d0 = d0.post,
                            d1 = d1.post,
                            d2 = d2.post,
                            d3 = d3.post,
                            d4 = d4.post),
              theta.h.post = theta.h.post,
             # eta.h.post = eta.h.post,
              Z.loc.post =Z.loc.post,
              Z.pred.post = Z.pred.post,
              lambda1.post = lambda1.post,
              lambda2.post = lambda2.post,
              lambda.post = lambda.post,
              mu.post = mu.post,
             cont = list(cont1, cont2, cont3, cont3b, cont4, cont4b))
         
  )
}



