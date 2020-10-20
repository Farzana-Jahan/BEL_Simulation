.libPaths("c:/software/Rpackages")

# irregular lattice data generation using sp package

# scottish lip cancer data
# data reading

# simulating underlying spatial random field
library(spdep)		# For generating random fields
library(ggplot2)
library(gridExtra)
library(rgdal)


map<- readOGR("C:/R dir/Leroux-Empirical-Likelihood-master/BEL_Simulation/Data/scotlip/scotlip.shp",verbose=F)

# Number of areas
N <- length(map)

# Compute centroids
Centroids <- matrix(NA, N, 2)
for(i in 1:N){
  bbox <- bbox(map@polygons[[i]])	# Bounding box of SLA
  Centroids[i,] <- c(bbox[1,1] + (bbox[1,2] - bbox[1,1]) / 2,
                     bbox[2,1] + (bbox[2,2] - bbox[2,1]) / 2)
}
Centroids <- data.frame(map$RECORD_ID[map$RECORD_ID], Centroids)
names(Centroids) <- c("Area_ID", "Long", "Lat")

# Matrix of distances
d <- matrix(NA, N, N)
Centroids <- Centroids[,2:3]
for(i in 1:N){
  d[i,] <- spDistsN1(as.matrix(Centroids), 
                     as.matrix(Centroids[i,]), longlat = FALSE)
}


n <- sqrt(N)        # For plots

# Decay functions  # is there any other way? # decay function: exponential # fields 
Gaus.decay <- function(d, str){
  exp(-0.5 * (d/str)^2)
}

# User-specific parameters (see Word document)
bandwidth <- 10		# Bandwith of decay function (Set by Master Script) # can have different values
r <- 3              # range of log-USRF to control effect size
y.total <- N*3      # Total observed counts  (Set by Master Script)  - I suggest base on N

#==========================================================================
  # Generate synthetic data # 
  #==========================================================================
Data_scot_sim<- list()
for(i in 1:10){
  set.seed(i)
  
  # 1) Generate log-underlying spatial random field
  log.USRF <- Gaus.decay(d, bandwidth)			# Spatial RF with specified autocorrelation
  log.USRF <- log.USRF %*% rnorm(N, 0.00, 0.01)			# Convert N x N SRF to N-length vector
  r.old <- max(log.USRF) - min(log.USRF)
  log.USRF <- log.USRF / r.old * r                # The range will now equal r
  # Although this doesn't need to be centred around zero, it shouuldn't be too far from zero
  if(max(log.USRF) > r*2/3){
    log.USRF <- log.USRF - max(log.USRF) + r/3
  }
  if(min(log.USRF) < -r*2/3){
    log.USRF <- log.USRF - min(log.USRF) - r/3
  }
  
  # 2) Generate observed values with covariate
  USRF <- exp(log.USRF)
  x<-  scale(runif(N)) # scaled covariate
  mean_y<-5+3*x # one covariate with intercept 5 and slope 3
  y<- rnorm(N,mean=mean_y)
# 3) Reorder y to match order of USRF
  log.RR <- log.USRF + mean_y  # mean_y calculated from sim covariate in step 2 
  RR <- exp(log.RR)
  ord.R <- order(RR)
  ord.y <- order(y)
  y <- y[ord.y][order(ord.R)]
  
  # without covariate 
  
  plot(density(y)) # y generated from a normal dstribution and then added some spatial effect/correlation
  
  # 3) Reorder y to match order of USRF
  ord.U <- order(USRF)
  ord.y <- order(y)
  y <- y[ord.y][order(ord.U)] # similar spatial pattern of y and USRF
  
  # 4) Compute expected values
  E <- y / exp(log.USRF)
  
  # 5) Rescale y and E
  # Problems to be solved:
  y <- y / sum(y) * y.total       # (y not discrete; sum(y) != sum(E))
  y <- round(y)                   # (sum(y) != sum(E))
  E <- E / sum(E) * sum(y)        # Done - all problems solved
  
  sum(y)
  sum(E)

  # created counts of y, underlying distribution is not actually Poisson
  
  #Combine data and export
  #==========================================================================
  
  Data_scot_sim[[i]] <- data.frame(
    Area.ID = 1:N,
    Observed = y,
    Expected = E,
    raw.SIR = y/E,
    x=x,
    Long = Centroids$Long,
    Lat = Centroids$Lat
  )
}
# saving 10 realisations of the data with moderate bandwith
write_csv(Data_scot_sim[[1]],"Data/sim_data/sim1_norm.csv")
write_csv(Data_scot_sim[[2]],"Data/scotlip/sim2_norm.csv")
write_csv(Data_scot_sim[[3]],"Data/scotlip/sim3_norm.csv")
write_csv(Data_scot_sim[[4]],"Data/scotlip/sim4_norm.csv")
write_csv(Data_scot_sim[[5]],"Data/scotlip/sim5_norm.csv")
write_csv(Data_scot_sim[[6]],"Data/scotlip/sim6_norm.csv")
write_csv(Data_scot_sim[[7]],"Data/scotlip/sim7_norm.csv")
write_csv(Data_scot_sim[[8]],"Data/scotlip/sim8_norm.csv")
write_csv(Data_scot_sim[[9]],"Data/scotlip/sim9_norm.csv")
write_csv(Data_scot_sim[[10]],"Data/scotlip/sim10_norm.csv")
save(Data_scot_sim,file="Data/sim_data/sim_norm_scotlip.RData")
# high spatial autocorrleation

