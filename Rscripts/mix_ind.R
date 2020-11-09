
# libraries
require(BELSpatial)
require(gmm)
require(emplik)
require(MASS)
require(ggplot2)        # For fortify(), ggplot()
require(readxl)         # For read_excel()
require(magrittr)       # For the pipe operator %>%
require(scales)         # For rescale()
require(dplyr)          # For inner_join(), bind_rows(),### between(), mutate()
require(gridExtra)      # For grid.arrange()
require(tidyr)          # For gather()
#require(spdep)




W<-readRDS("Data/sim_data/W_25areas")

ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

load("Data/sim_data/sim_mix_25.RData")


BEL_ind_mix_25<-list()

for(i in 1:5){
  data1<-Data_mix[[i]]
  data1$raw.SIR[data1$raw.SIR==0]<-0.1
  y<- log(data1$raw.SIR)
  x<- cbind(1, data1$x)
  # initial values needed before fitting models
  
  n<- length(data1$raw.SIR) # no. of observations
  p<- 2 # no. of covariates
  alpha_1<-1 # hyperparamter for tau prior
  alpha_2<-0.01 # hyperparamter for tau prior
  tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
  tau_init<- 1/tau_inv_init
  g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
  prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
  beta_init<- rnorm(2,prior_mean_beta, (1/g)*tau_inv_init)
  wi_init<- 1/n # y be the response variable from the data
  psi_init <- rep(0,n)
  
  
  # calculating MELE of Beta, beta_mele
  var<- as.numeric(var(y- x%*%beta_init))
  wi=wi_init
  beta_mele<- mele(x,tet=beta_init,y=y,var=var)
  mu_init<- x%*% beta_mele + psi_init
  beta_init<-beta_mele
  wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
  wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
  wi<-wi_mu
  
  # fitting BEL BYM model taking rho= 1
  library(parallel)
  cluster<-makeCluster(3)
  #clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
  clusterEvalQ(cl=cluster,library(BELSpatial))
  clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))


BEL_ind_mix_25[[i]]<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0,niter=1000000,
                                                                              beta_init, psi_init, tau_init,R, wi, sd_psi=1.8, 
                                                                              sd_beta=2, sd_tau=0.6)})
}
save(BEL_ind_mix_25,file="Results/BEL_BYM_mix_25.RData")
