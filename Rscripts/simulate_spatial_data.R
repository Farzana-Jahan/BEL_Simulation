.libPaths("c:/software/Rpackages")

# regular lattice data generation using sp package



# simulating underlying spatial random field
library(spdep)		# For generating random fields
library(ggplot2)
library(gridExtra)
library(rgdal)


map<- readOGR("Data/sim_data/My New Shapefile_100areas.shp",verbose=F)

# Number of areas
N <- length(map)

# Compute centroids
Centroids <- matrix(NA, N, 2)
for(i in 1:N){
  bbox <- bbox(map@polygons[[i]])	# Bounding box of SLA
  Centroids[i,] <- c(bbox[1,1] + (bbox[1,2] - bbox[1,1]) / 2,
                     bbox[2,1] + (bbox[2,2] - bbox[2,1]) / 2)
}
Centroids <- data.frame(map$id[map$id], Centroids)
names(Centroids) <- c("Area_ID", "Long", "Lat")

# Matrix of distances
d <- matrix(NA, N, N)
Centroids <- Centroids[,2:3]
for(i in 1:N){
  d[i,] <- spDistsN1(as.matrix(Centroids), 
                     as.matrix(Centroids[i,]), longlat = FALSE)
}


n <- sqrt(N)        # For plots

# Decay functions  
Gaus.decay <- function(d, str){
  exp(-0.5 * (d/str)^2)
}

# User-specific parameters (see Word document)
bandwidth <- 10		# Bandwith of decay function (Set by Master Script) # can have different values
r <- 3              # range of log-USRF to control effect size
y.total <- N*3      # Total observed counts  (Set by Master Script)  - I suggest base on N

#==========================================================================
  # Generate synthetic data # # scenario 1, area=100, N small, count (smaller),bandwith high 10,Normal counts
  #==========================================================================
Data_sim1<- list()
for(i in 1:5){
  set.seed(i)
  
  # 1) Generate log-underlying spatial random field
  log.USRF <- Gaus.decay(d, bandwidth)			# Spatial RF with specified autocorrelation
  log.USRF <- log.USRF %*% rnorm(N, 0.00, 0.01)			# Convert N x N SRF to N-length vector
  r.old <- max(log.USRF) - min(log.USRF)
  log.USRF <- log.USRF / r.old * r                # The range will now equal r
  # Although this doesn't need to be centred around zero, it shouuldn't be too far from zero
  if(max(log.USRF) > r*(2/3)){
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
  
  Data_sim1[[i]] <- data.frame(
    Area.ID = 1:N,
    Observed = y,
    Expected = E,
    raw.SIR = y/E,
    x=x,
    Long = Centroids$Long,
    Lat = Centroids$Lat
  )
}


save(Data_sim1,file="Data/sim_data/sim_norm_100.RData",version = 2)
# high spatial autocorrleation, bandwith=5, 100 areas, norm
#==========================================================================
# Generate synthetic data # # scenario 1, area=100, N small, count (smaller),bandwith moderate,Normal counts
#==========================================================================
Data_sim2<- list()
for(i in 1:5){
  set.seed(i)
  
  # 1) Generate log-underlying spatial random field
  log.USRF <- Gaus.decay(d, 5)			# Spatial RF with specified autocorrelation
  log.USRF <- log.USRF %*% rnorm(N, 0.00, 0.01)			# Convert N x N SRF to N-length vector
  r.old <- max(log.USRF) - min(log.USRF)
  log.USRF <- log.USRF / r.old * r                # The range will now equal r
  # Although this doesn't need to be centred around zero, it shouuldn't be too far from zero
  if(max(log.USRF) > r*(2/3)){
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
  
  Data_sim2[[i]] <- data.frame(
    Area.ID = 1:N,
    Observed = y,
    Expected = E,
    raw.SIR = y/E,
    x=x,
    Long = Centroids$Long,
    Lat = Centroids$Lat
  )
}


save(Data_sim2,file="Data/sim_data/sim2_norm_100.RData",version = 2)
#==========================================================================
# Generate synthetic data # # scenario 1, area=100, N small, count (smaller),bandwith very small,Normal counts
#==========================================================================
Data_sim3<- list()
for(i in 1:5){
  set.seed(i)
  
  # 1) Generate log-underlying spatial random field
  log.USRF <- Gaus.decay(d, 1)			# Spatial RF with specified autocorrelation
  log.USRF <- log.USRF %*% rnorm(N, 0.00, 0.01)			# Convert N x N SRF to N-length vector
  r.old <- max(log.USRF) - min(log.USRF)
  log.USRF <- log.USRF / r.old * r                # The range will now equal r
  # Although this doesn't need to be centred around zero, it shouuldn't be too far from zero
  if(max(log.USRF) > r*(2/3)){
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
  E <- round(y / exp(log.USRF),0)
  
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
  
  Data_sim3[[i]] <- data.frame(
    Area.ID = 1:N,
    Observed = y,
    Expected = E,
    raw.SIR = y/E,
    x=x,
    Long = Centroids$Long,
    Lat = Centroids$Lat
  )
}


save(Data_sim3,file="Data/sim_data/sim3_norm_100.RData",version = 2)

