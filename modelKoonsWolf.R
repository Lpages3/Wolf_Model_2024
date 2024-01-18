# Supplement for "A niche for null models in 
# Adaptive Resource Management"

###########################################################
# We fit the null model of persistence with a random walk 
# to available midcontinent mallard abundance data using
# a Bayesian state-space model. Upon fitting the model
# recursively each year, we made 1-step-ahead forecasts 
# that could later be compared to data collected in 
# the subsequent year.
###########################################################
sink("persistnull.jags")
cat("
model {

#--------------------------------------------------------
# 1. Define the priors for the parameters
#--------------------------------------------------------
# Informative prior for initial abundance based on the 
# data and associated s.e.
N[1] ~ dnorm(17, pow(0.28,-2))I(0,) 

# Prior for the temporal process variance
# note that tauPro is 1/variance
sigPro ~ dunif(0,5)
tauPro <- pow(sigPro,-2)
lambda ~ dunif(0, 5)

#--------------------------------------------------------
# 2. Process model of state-space likelihood
#--------------------------------------------------------
for (t in 1:years) {
  N[t+1] <- N[t]*lambda
}

#--------------------------------------------------------
# 3. Observation model of state-space likelihood
#--------------------------------------------------------
for (t in 1:years) {
  y[t] ~ dnorm(N[t], pow(yse[t],-2)) # both y and yse are
                                     # provided as data
}

}
",fill = TRUE)
sink()

library(tidyverse)

harvest <- c(0,0,0,0,0,1,0,0,2,1,2,0,0,1,0,4,4,6,18,36,34,42,51,98,105,103,169)
CMR = c(17.1,35.4, 47.7,25.1,62.6,47.9,81.7,110.5,102.7,135.9,132.6,101.7,130.3,141.4,141.5,175.5,210.3,174.5,353.6,280.2,376.7,561.2,571.9,682.4,645.7,783.8,868)
thedata <- cbind(round(CMR), harvest)
colnames(thedata) <- c("N", "H")
thedata <- as.data.frame(thedata)

Y=thedata$N
Yse[26]=0.3
Yse[27]=0.28
# Description of data:
# Y: annual estimates of midcontinent mallard abundance (in millions)
# Yse: s.e. of Y

# set up iterative forecasting with each 
# successive year of monitoring
library(jagsUI)
perspredict <- matrix(NA,length(Y),1)
perspredict[1] <- Y[1]
for (i in 2:length(Y)){
  # Bundle data
  y = Y[1:i]
  yse = Yse[1:i]
  years = i 
  jags.data <- list(y = y, yse = yse, years = years)
  
  # Initial values
  inits <- function(){list(sigPro = runif(1,0.7,1),
                           N = c(runif(1,7,8.5), rep(NA, years)) )}
  
  # Parameters monitored
  parameters <- c("sigPro","N")
  
  # MCMC settings
  ni <- 50000; nt <- 15; na <- 5000; nb <- 20000; nc <- 3
  
  # Call JAGS from R
  persistnull <- jags(jags.data, inits, parameters, "persistnull.jags", 
                      n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.adapt = na, n.burnin = nb)
  
  perspredict[i] <- median(persistnull$sims.list$N[,(i+1)])
}

ggplot()+
  geom_line(aes(x = 1995 + 1:unique(years), y = perspredict), lty="dashed", color="red") +
  geom_point(aes(x = 1995 + 1:unique(years), y = Y), color="black")+
  labs(x = "Years", y = "Annual estimates of midcontinent mallard abundance (in millions)")
