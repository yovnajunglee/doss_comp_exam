
## Code adapted from Fonseca et al. (2023) 
## https://github.com/thaiscofonseca/DynGLG



#----------------------------------
### useful packages
require(MCMCpack)
require(mvtnorm)
require("MASS")
require(msm)
library(coda)


#----------------------------------
### distance in space

distf = function(coord){
  nI = dim(coord)[1]
  result = diag(nI)*0
  for (i in 1:nI){ for (j in i:nI) {
    result[i,j] = sqrt((coord[i,1]-coord[j,1])^2+(coord[i,2]-coord[j,2])^2)
    result[j,i] = result[i,j]
  }}
  result
}


#-------------------------------------
#### correlation matrix

## using cauchy correlation 
CauchyCov<- function(distance,phi,alpha){
  value1 = (distance*(1/phi)) ^ alpha;
  (1+value1) ^ (-1)
}


#------------------------------------------------------------------------
#------------------------------------------------------------------------

######################################################################### FFBS
# Forward Filtering Backward Sampling for theta - mean
FFBS= function(z,V,discount,m0,C0,Ftt,Gt, loglamb1, loglamb2, sigma2, tau2, In, h = 10){
  
  r <- ncol(z)   ## 66 stations
  T <- nrow(z)   ## 31 times
  n <- nrow(Ftt[[1]])   ## 4 (we have 3 covariates plus 1 vector)
  
  # arrays
  mt <- array(data = NA,dim = c(n,1,T)) 
  Ct <- array(data = NA,dim = c(n,n,T))  
  at <- array(data = NA,dim = c(n,1,T)) 
  Rt <- array(data = NA,dim = c(n,n,T)) 
  
  
  ft <- array(data = NA,dim = c(r,1,T))   
  Qt <- array(data = NA,dim = c(r,r,T))  
  At <- array(data = NA,dim = c(n,r,T))  
  et <- array(data = NA,dim = c(r,1,T)) 
  
  
  B  = array(rep(diag(n),T),dim=c(n,n,T))
  H  = array(rep(diag(n),T),dim=c(n,n,T))
  
  m.t =  array(data = NA,dim = c(n,1,T))
  C.t = array(rep(diag(n),T),dim=c(n,n,T))
  delta  = array(data = NA,dim = c(n,1,T))
  
  
  #### run Kalman filter - petris-petroni page
  
  # Forward filtering
  for(t in 1:T){
    ## b
    #[theta_t|D_t-1]
    if(t == 1){
      W<- C0*(1 - discount) / discount
      at[,,t] <- m0
      Rt[,,t] <- Gt%*%C0%*%t(Gt) + W
    }else{
      W<- Ct[,,t-1]*(1- discount)/discount
      at[,,t] <- Gt%*%mt[,,t-1]
      Rt[,,t] <- Gt%*%Ct[,,t-1]%*%t(Gt) + W
    }
    
    ## c
    #[y_t|D_t-1]
    Ft = Ftt[[t]]
    ft[,,t] <- t(Ft)%*%at[,,t]
    lamb1 = as.vector(exp(loglamb1)) 
    lamb2 = exp(loglamb2[t])
    Sig2.star = sigma2*V #calc. observational variance
    Sig2 = (diag(1/sqrt(lamb1*lamb2))%*%Sig2.star%*%diag(1/sqrt(lamb1*lamb2)))
    Qt[,,t] <- t(Ft)%*%Rt[,,t]%*%(Ft) + Sig2
    ## d
    # [theta_t|D_t]
    At[,,t] <- Rt[,,t]%*%Ft%*%solve(Qt[,,t])
    et[,,t] <- z[t,] - ft[,,t]
    
    mt[,,t] <- at[,,t] + At[,,t]%*%et[,,t]   ## posterior mean
    Ct[,,t] <- Rt[,,t] - At[,,t]%*%(Qt[,,t])%*%t(At[,,t])  ### posterior variance
    
  }
  
  #go backwards
  m.t[,,T] = mt[,,T]
  C.t[,,T] = Ct[,,T]
  
  #step 2 pg 135
  delta[,,T]  = mvrnorm(1, m.t[,,T], C.t[,,T])   ## draw deltaT ~ N(mt, Ct)
  
  ### petris-petroni page 162  (step 3)
  
  # Backward smoothing
  # [theta_t|D_T]
  #Gt fixed for all t
  ### petris-petroni page 162  (step 3) [theta_t|theta_t+1, D_t]
  
  for(t in (T-1):1){
    B[,,t] = Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])
    
    H[,,t] = Ct[,,t] - B[,,t]%*%Rt[,,t+1]%*%t(B[,,t])
    
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.10 in book but put + before B
    # notation m.t = a_T(t-T)
    m.t[,,t] <- mt[,,t] + Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])%*%(m.t[,,t+1] - at[,,t+1])
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.11 in book but put + before B
    # notation C.t = R_T(t-T)
    C.t[,,t] <- Ct[,,t] - Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - C.t[,,t+1])%*%t(Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1]))
    
    
    # do we actually need m.t, C.t??
    
    delta[,,t]  = mvrnorm(1,mt[,,t]+B[,,t]%*%(delta[,,t+1]-at[,,t+1]),H[,,t])   ##   ht = mt[,,t]+B[,,t]%*%(delta[,,t+1]-at[,,t+1])
  }
  
  
  # Predictive distribution of theta_t+h
  # pg 124
  # eq. 4.9
  a.h = array(data = NA,dim = c(n,1,h))
  R.h = array(data = NA,dim = c(n,n,h))
  theta.h  = array(data = NA,dim = c(n,1,h))
  
  WT<- Ct[,,T]*(1- discount)/discount
  
  for(hh in 1:h){
    
    if(hh==1){
      a.h[,,hh] = Gt%*%mt[,,T]
      R.h[,,hh] = Gt%*%Ct[,,T]%*%t(Gt) + WT
    }else{
      a.h[,,hh] = Gt%*%a.h[,,hh-1]
      R.h[,,hh] = Gt%*%R.h[,,hh-1]%*%t(Gt) + WT
    }
    theta.h[,,hh] = mvrnorm(1,a.h[,,hh],R.h[,,hh])
    
  }
  
  return(list(d=delta,m=mt,C=Ct,m.m=m.t,C.C=C.t,theta.h=theta.h))
}


FFBS.g= function(z,V,discount,m0,C0,Ftt,Gt, loglamb1, loglamb2, sigma2, tau2, In, h = 10){
  
  r <- ncol(z)   ## 66 stations
  T <- nrow(z)   ## 31 times
  n <- nrow(Ftt[[1]])   ## 4 (we have 3 covariates plus 1 vector)
  
  # arrays
  mt <- array(data = NA,dim = c(n,1,T)) 
  Ct <- array(data = NA,dim = c(n,n,T))  
  at <- array(data = NA,dim = c(n,1,T)) 
  Rt <- array(data = NA,dim = c(n,n,T)) 
  
  
  ft <- array(data = NA,dim = c(r,1,T))   
  Qt <- array(data = NA,dim = c(r,r,T))  
  At <- array(data = NA,dim = c(n,r,T))  
  et <- array(data = NA,dim = c(r,1,T)) 
  
  
  B  = array(rep(diag(n),T),dim=c(n,n,T))
  H  = array(rep(diag(n),T),dim=c(n,n,T))
  
  m.t =  array(data = NA,dim = c(n,1,T))
  C.t = array(rep(diag(n),T),dim=c(n,n,T))
  delta  = array(data = NA,dim = c(n,1,T))
  
  
  #### run Kalman filter - petris-petroni page
  
  # Forward filtering
  for(t in 1:T){
    ## b
    #[theta_t|D_t-1]
    if(t == 1){
      W<- C0*(1 - discount) / discount
      at[,,t] <- m0
      Rt[,,t] <- Gt%*%C0%*%t(Gt) + W
    }else{
      W<- Ct[,,t-1]*(1- discount)/discount
      at[,,t] <- Gt%*%mt[,,t-1]
      Rt[,,t] <- Gt%*%Ct[,,t-1]%*%t(Gt) + W
    }
    
    ## c
    #[y_t|D_t-1]
    Ft = Ftt[[t]]
    ft[,,t] <- t(Ft)%*%at[,,t]
    lamb1 = as.vector(exp(loglamb1)) 
    lamb2 = exp(loglamb2[t])
    Sig2.star = sigma2*V #calc. observational variance
    Sig2 = Sig2.star
    #(diag(1/sqrt(lamb1*lamb2))%*%Sig2.star%*%diag(1/sqrt(lamb1*lamb2)))
    Qt[,,t] <- t(Ft)%*%Rt[,,t]%*%(Ft) + Sig2
    ## d
    # [theta_t|D_t]
    At[,,t] <- Rt[,,t]%*%Ft%*%solve(Qt[,,t])
    et[,,t] <- z[t,] - ft[,,t]
    
    mt[,,t] <- at[,,t] + At[,,t]%*%et[,,t]   ## posterior mean
    Ct[,,t] <- Rt[,,t] - At[,,t]%*%(Qt[,,t])%*%t(At[,,t])  ### posterior variance
    # Symmetrise the matrix
    Ct[,,t]=0.5*(Ct[,,t]+t(Ct[,,t]))
    #print(paste(t,isSymmetric(Ct[,,t])))
    
  }
  
  #go backwards
  m.t[,,T] = mt[,,T]
  C.t[,,T] = Ct[,,T]
  
  #step 2 pg 135
  delta[,,T]  = mvrnorm(1, m.t[,,T], C.t[,,T])   ## draw deltaT ~ N(mt, Ct)
  
  ### petris-petroni page 162  (step 3)
  
  # Backward smoothing
  # [theta_t|D_T]
  #Gt fixed for all t
  ### petris-petroni page 162  (step 3) [theta_t|theta_t+1, D_t]
  
  for(t in (T-1):1){
    B[,,t] = Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])
    
    H[,,t] = Ct[,,t] - B[,,t]%*%Rt[,,t+1]%*%t(B[,,t])
    
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.10 in book but put + before B
    # notation m.t = a_T(t-T)
    m.t[,,t] <- mt[,,t] + Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])%*%(m.t[,,t+1] - at[,,t+1])
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.11 in book but put + before B
    # notation C.t = R_T(t-T)
    C.t[,,t] <- Ct[,,t] - Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - C.t[,,t+1])%*%t(Ct[,,t]%*%t(Gt)%*%solve(Rt[,,t+1]))
    
    
    # do we actually need m.t, C.t??
    delta[,,t]  = mvrnorm(1,mt[,,t]+B[,,t]%*%(delta[,,t+1]-at[,,t+1]),H[,,t])   ##   ht = mt[,,t]+B[,,t]%*%(delta[,,t+1]-at[,,t+1])
  }
  
  
  # Predictive distribution of theta_t+h
  # pg 124
  # eq. 4.9
  a.h = array(data = NA,dim = c(n,1,h))
  R.h = array(data = NA,dim = c(n,n,h))
  theta.h  = array(data = NA,dim = c(n,1,h))
  
  WT<- Ct[,,T]*(1- discount)/discount
  
  for(hh in 1:h){
    
    if(hh==1){
      a.h[,,hh] = Gt%*%mt[,,T]
      R.h[,,hh] = Gt%*%Ct[,,T]%*%t(Gt) + WT
    }else{
      a.h[,,hh] = Gt%*%a.h[,,hh-1]
      R.h[,,hh] = Gt%*%R.h[,,hh-1]%*%t(Gt) + WT
    }
    theta.h[,,hh] = mvrnorm(1,a.h[,,hh],R.h[,,hh])
    
  }
  
  return(list(d=delta,m=mt,C=Ct,m.m=m.t,C.C=C.t,theta.h=theta.h))
}

############################################################
# Forward Filtering Backward Sampling for mu
FFBS.mu= function(l.lamb,V1,d2,m1,C11,Ftt,Gt1, h = 10){
  
  r <- ncol(l.lamb)   ## 66 stations
  T <- nrow(l.lamb)   ## 31 times
  #  n <- nrow(Ftt[[1]])   ## 4 (we have 3 covariates plus 1 vector)
  n <- nrow(Ftt)   
  
  # arrays
  mt <- array(data = NA,dim = c(n,1,T)) 
  Ct <- array(data = NA,dim = c(n,n,T))  
  at <- array(data = NA,dim = c(n,1,T)) 
  Rt <- array(data = NA,dim = c(n,n,T)) 
  ft <- array(data = NA,dim = c(r,1,T))   
  Qt <- array(data = NA,dim = c(r,r,T))  
  At <- array(data = NA,dim = c(n,r,T))  
  et <- array(data = NA,dim = c(r,1,T)) 
  
  
  B  = array(rep(diag(n),T),dim=c(n,n,T))
  H  = array(rep(diag(n),T),dim=c(n,n,T))
  
  m.t =  array(data = NA,dim = c(n,1,T))
  C.t = array(rep(diag(n),T),dim=c(n,n,T))
  delta  = array(data = NA,dim = c(n,1,T))
  
  
  #### run Kalman filter - petris-petroni / west-harrison
  for(t in 1:T){
    ## b
    if(t == 1){
      W1<- C11*(1 - d2) / d2
      at[,,t] <- m1
      Rt[,,t] <- Gt1%*%C11%*%t(Gt1) + W1
    }else{
      W1<- Ct[,,t-1]*(1- d2)/d2
      at[,,t] <- Gt1%*%mt[,,t-1]
      Rt[,,t] <- Gt1%*%Ct[,,t-1]%*%t(Gt1) + W1
    }
    
    ## c
    Ft1 = Ftt[,t]
    ft[,,t] <- t(Ft1)%*%at[,,t]
    Qt[,,t] <- t(Ft1)%*%Rt[,,t]%*%(Ft1) + V1
    
    ## d
    At[,,t] <- Rt[,,t]%*%Ft1%*%solve(Qt[,,t])
    et[,,t] <- l.lamb[t,] - ft[,,t]
    
    mt[,,t] <- at[,,t] + At[,,t]*et[,,t]   ## posterior mean
    Ct[,,t] <- Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[,,t]  ### posterior variance
    
    
  }
  
  m.t[,,T] = mt[,,T]
  C.t[,,T] = Ct[,,T]
  delta[,,T]  = mvrnorm(1, m.t[,,T],C.t[,,T])   ## draw deltaT ~ N(mt, Ct)
  
  ### petris-petroni page 162  (step 3) [theta_t|theta_t+1, D_t]
  for(t in (T-1):1){
    B[,,t] = Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1])
    
    H[,,t] = Ct[,,t] - B[,,t]%*%Rt[,,t+1]%*%t(B[,,t]) # eq. 4.13
    
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.10 in book but put + before B
    # notation m.t = a_T(t-T)
    m.t[,,t] <- mt[,,t] + Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1])%*%(m.t[,,t+1] - at[,,t+1])
    
    # couldve replaced Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]) by  B[,,t] 
    # same as eq. 4.11 in book but put + before B
    # notation C.t = R_T(t-T)
    C.t[,,t] <- Ct[,,t] - Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - C.t[,,t+1])%*%t(Ct[,,t]%*%t(Gt1)%*%solve(Rt[,,t+1]))
    
    # do we actually need m.t, C.t??
    delta[,,t]  = mvrnorm(1,mt[,,t]+B[,,t]%*%(delta[,,t+1]-at[,,t+1]), H[,,t])  
    
  }
  
  # Predictive distribution of theta_t+h
  # pg 124
  # eq. 4.9
  a.h = array(data = NA,dim = c(n,1,h))
  R.h = array(data = NA,dim = c(n,n,h))
  eta.h  = array(data = NA,dim = c(n,1,h))
  
  WT<- Ct[,,T]*(1- d2)/d2
  
  for(hh in 1:h){
    
    if(hh==1){
      a.h[,,hh] = Gt1%*%mt[,,T]
      R.h[,,hh] = Gt1%*%Ct[,,T]%*%t(Gt1) + WT
    }else{
      a.h[,,hh] = Gt1%*%a.h[,,hh-1]
      R.h[,,hh] = Gt1%*%R.h[,,hh-1]%*%t(Gt1) + WT
    }
    eta.h[,,hh] = mvrnorm(1,a.h[,,hh],R.h[,,hh])
  }
  
  return(list(d=delta,m=mt,C=Ct,m.m=m.t,C.C=C.t, eta.h = eta.h))
}

#************************** log-likelihood *******************
#Log likelihood of eq. (1)
LogVeroNG = function(vecz,distanceS,Plat,Plong,Ptemp,d0,d1,d2,d3,d4,sig2,phi,alpha,loglamb1,loglamb2, tau2){
  I = length(loglamb1)
  J = length(loglamb2)
  lvero=NULL
  det<- NULL
  lamb1 = exp(loglamb1)
  C1 = CauchyCov(distanceS,phi,alpha)
  # Precision matrix ( without lambda 2)
  Cinv = diag(as.vector(sqrt(lamb1)))%*%solve(sig2*C1)%*%diag(as.vector(sqrt(lamb1)))
  # ==  Determinant of Cinv
  Cdet = (determinant(diag(as.vector(1/sqrt(lamb1)))%*%(sig2*C1)%*%diag(as.vector(1/sqrt(lamb1))),logarithm = TRUE)$modulus)[[1]]
  for(t in 1:J){
    lamb2t = exp(loglamb2[t])
    Sig2inv = lamb2t*Cinv
    mu<- d0[t] +  d1[t]*Plat +  d2[t]*Plong +  d3[t]*Ptemp[t,]
    lvero[t] =  t(vecz[t,]-mu)%*%Sig2inv%*%(vecz[t,]-mu)
    det[t]<- -I*log(lamb2t)+Cdet
  }
  l.vero<-    -0.5*sum(det)- 0.5*sum(lvero)
  l.vero[[1]]
}


#Log likelihood of eq. (1) for region
LogVeroNGt = function(vecz,distanceS,Plat,Plong,Ptemp,d0,d1,d2,d3,d4,sig2,phi,alpha,loglamb1,loglamb2, tau2,region){
  
  I = length(loglamb1)
  J = length(loglamb2)
  lvero=rep(NA,J)
  det = rep(NA,J)
  I.aux = diag(1, length(vecz[1,]))
  lamb1 = exp(loglamb1)
  C1 = CauchyCov(distanceS,phi,alpha)
  Cinv = diag(as.vector(sqrt(lamb1)))%*%solve(sig2*C1)%*%diag(as.vector(sqrt(lamb1)))
  Cdet = (determinant(diag(as.vector(1/sqrt(lamb1)))%*%(sig2*C1)%*%diag(as.vector(1/sqrt(lamb1))),logarithm = TRUE)$modulus)[[1]]
  for(t in region){
    lamb2t = exp(loglamb2[t])
    Sig2inv = lamb2t*Cinv
    mu<- d0[t] +  d1[t]*Plat +  d2[t]*Plong +  d3[t]*Ptemp[t,]
    lvero[t] =  t(vecz[t,]-mu)%*%Sig2inv%*%(vecz[t,]-mu)
    det[t]<- -I*log(lamb2t)+Cdet
  }
  l.vero<-    -0.5*sum(det[region])- 0.5*sum(lvero[region])
  l.vero[[1]]
}


LogVeroLogLamb1 = function(loglamb1,distanceS,phi,alpha,nu1){ 
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = nu1*exp(-distanceS/phi)  ## covariance matrix for lambda
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}

LogVeroLogLamb2 = function(loglamb2,nu2,mut){
  J = length(loglamb2)
  lverolamb<- NULL
  Wlamb = nu2*1   ## covariance matrix for lambda
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu2/2*rep(1, J)+mut
  lverolamb =  sum(dnorm(loglamb2,mu.loglamb,sqrt(Wlamb),log=T))
  as.vector(lverolamb)
}



# =====



####################### here read functions

# Log-likelihood for log(lambda_1)
# Regression for spatial mixing process
LogVeroLogLamb1reg = function(loglamb1,distanceS,phi,alpha,nu1,beta1,x1){ 
  # Eq. (3) in paper
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = nu1*exp(-distanceS/phi)  ## covariance matrix for lambda1
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)+x1%*%t(beta1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}
######### main code

####################### here read functions

# Log-likelihood for log(lambda_1)
# Regression for spatial mixing process
LogVeroLogLamb1maternreg = function(loglamb1,distanceS,kappa,nu,nu1,beta1,x1){ 
  # Eq. (3) in paper
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = matern.covariance(distanceS, kappa, nu, nu1)
  #exp(-distanceS/phi)  ## covariance matrix for lambda1
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)+x1%*%t(beta1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}

LogVeroLogLamb1maternregfields = function(loglamb1,distanceS,kappa,nu,nu1,beta1,x1){ 
  # Eq. (3) in paper
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = Matern(distanceS, range = kappa, smoothness = nu, phi = nu1)
  #exp(-distanceS/phi)  ## covariance matrix for lambda1
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)+x1%*%t(beta1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}
######### main code

LogVeroLogLamb1matern = function(loglamb1,distanceS,kappa,nu,nu1){ 
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = matern.covariance(distanceS, kappa, nu, nu1)
  #nu1*exp(-distanceS/phi)  ## covariance matrix for lambda
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}



LogVeroLogLamb1maternfields = function(loglamb1,distanceS,kappa,nu,nu1){ 
  I1 = length(loglamb1)
  lverolamb<- NULL
  Wlamb = Matern(distanceS, range = kappa, smoothness = nu, phi = nu1)
  #nu1*exp(-distanceS/phi)  ## covariance matrix for lambda
  W.aux<- solve(Wlamb)
  mu.loglamb = -nu1/2*rep(1, I1)
  lverolamb =  t(loglamb1-mu.loglamb)%*%W.aux%*%(loglamb1-mu.loglamb)
  l.vero<-    -0.5*determinant(Wlamb,logarithm = TRUE)$modulus- 0.5*sum(lverolamb)  
  l.vero[[1]]
  
}
