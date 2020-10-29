# for office PC only run the /libpaths code, otherwise comment it out
# .libPaths("c:/software/Rpackages")
#devtools::install_github("danwkenn/BELSpatial") # installing the package
#devtools::install_local("C:\\R dir\\Leroux-Empirical-Likelihood-master\\Leroux-Empirical-Likelihood-master\\BELSpatial",force=TRUE)
# loading the data (Scottish Lip Cancer Data)

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

#W<- read.gal("Data/sim_data/My New Shapefile_2500areas.gal")
#W<-nb2mat(W,style="B")
#creating symmetric neighbourhood matrix for BYM in CARBAYES
#rownames(W)<-c()
#ind <- upper.tri(W)
#W[ind] <- t(W)[ind] 
#isSymmetric(W)


W<-readRDS("Data/sim_data/W_2500areas.RDS")

ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

load("Data/sim_data/sim_norm_2500.RData")
# creating blank list for all the diff realisations of the simulations
BEL_BYM_sim_2500<- list()
BEL_Leroux_sim_2500<-list()
Porter_BSHEL_sim_2500<-list()
BEL_ind_sim_2500<-list()

for(i in 1:5){
  data1<-Data_sim1[[i]]
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
  
  BEL_ind_sim_2500[[i]]<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0,niter=25000000,
                                                              beta_init, psi_init, tau_init,R, wi, sd_psi=1.6, 
                                                                  sd_beta=3, sd_tau=0.3)})
  
  
  # fitting BEL BYM model taking rho= 0.95
  clusterEvalQ(cl=cluster,library(BELSpatial))
  clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
  
  BEL_BYM_sim_2500<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=1, niter=25000000,
                                                                                     beta_init, psi_init, tau_init,R, wi, sd_psi=0.2, 
                                                                                     sd_beta=1, sd_tau=0.3)})
  
  clusterEvalQ(cl=cluster,library(BELSpatial))
  clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
  
  BEL_Leroux_sim_2500<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.95, niter=25000000,
                                                                              beta_init, psi_init, tau_init,R, wi, sd_psi=0.2, 
                                                                              sd_beta=1.3, sd_tau=0.3)})
  
  
  # Porter's BSHEL model
  B<-W
  B_plus<-diag(rowSums(B))
  M=M_create(y,x,B)
  MBM=MBM_create(M,B,B_plus)
  q=dim(MBM)[2]
  psi_init <- rep(0,q) 
  wi=wi_init
  var<- as.numeric(var(y- x%*%beta_init))
  beta_mele<- mele( x = x, tet= beta_init,y=y,var=var) 
  mu_init<- x%*% beta_mele + M%*%psi_init
  beta_init<-beta_mele
  
  wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
  wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
  wi<-wi_mu
  
  #Fitting the Porter BSHEL model
  
  #clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
  clusterEvalQ(cl=cluster,library(BELSpatial))
  clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init"
                                       ,"B","B_plus","q","M","MBM", "wi"))
  Porter_BSHEL_sim_2500[[i]]<-clusterApply(cl=cluster, x=1:3, fun= function(z){BSHEL(y,x,n,p,q,var,niter=25000000,beta_init, 
                                                                                   psi_init, tau_init,M,MBM, wi, 
                                                                                   sd_psi=0.0000001, 
                                                                                   sd_beta=0.00000001, sd_tau=1.5)})
  
  
}

save(BEL_BYM_sim_2500,file="Results/BEL_BYM_sim_2500_sim1.RData")
save(BEL_ind_sim_2500,file="Results/BEL_BYM_sim_2500_sim1.RData")
save(BEL_leroux_sim_2500,file="Results/BEL_Leroux_sim_2500_sim1.RData")
save(Porter_BSHEL_sim_2500,file="Results/BSHEL_porter_sim_2500_sim1.RData")

