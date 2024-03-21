##--- Finite-mixture capture-recapture models
##--- O. Gimenez & D. Turek, December 2019 & October 2020

##--- Check out our paper
##--- Turek, D., C. Wehrhahn, Gimenez O. (2021). Bayesian Non-Parametric Detection Heterogeneity in Ecological Models. 
##--- Environmental and Ecological Statistics 28: 355-381.
##--- PDF available here: https://arxiv.org/abs/2007.10163

##--- More on heterogeneity in capture-recapture models in
##--- https://onlinelibrary.wiley.com/doi/abs/10.1111/oik.04532
##--- PDF available here: https://bit.ly/35KfCyj

## load packages
# library(nimble)
# library(MCMCvis)

library(R2jags)

##--- 1. simulate data

## simulate detection for each individual, 
## by first assigning them to the class 1 or 2, 
## then using the corresponding detection 

p_class1 <- 0.7 # detection class 1
p_class2 <- 0.4 # detection class 2
prop_class1 <- 0.6 # pi
phi <- 0.7 # survival 
nind <- 200 # nb of ind
first <- rep(1,nind) # date of first capture
nyear <- 5 # duration of the study
z <- data <- y <- matrix(NA,nrow=nind,ncol=nyear)
p <- rep(NA,nind)

which_mixture <- rep(NA,nind)
for(i in 1:nind) {
    which_mixture[i] <- rbinom(1,1,prop_class1) # assign ind i to a class with prob pi
    if(which_mixture[i] == 1) {
        p[i] <- p_class1
    } else { 
        p[i] <- p_class2
    }
}

## simulate the encounter histories:
for(i in 1:nind) {
    z[i,first[i]] <- y[i,first[i]] <- 1
    for(j in (first[i]+1):nyear) {
        z[i,j] <- rbinom(1,1,phi*z[i,j-1])
        y[i,j] <- rbinom(1,1,z[i,j]*p[i])
    }
}
y[is.na(y)] <- 0

## generate first:
get.first <- function(x) min(which(x!=0)) # get occasion of marking 
first <- apply(y, 1, get.first) 


##--- 2. Capture-recapture model with constant parameters 
##---    and heterogeneity in detection using 2-class finite mixtures

## Model fitting

code1 <- function(){
    for(i in 1:nind) { # for each ind
        group[i] ~ dbern(pi)
        z[i,first[i]] <- 1   
        y[i,first[i]] ~ dbern(z[i,first[i]]) 
        ##            
        for(j in (first[i]+1):nyear) { # loop over time
            ## STATE EQUATIONS ##
            ## draw states at j given states at j-1
            z[i,j] ~ dbern(phi * z[i,j-1])  
            ## OBSERVATION EQUATIONS ##
            ## draw observations at j given states at j
            y[i,j] ~ dbern(z[i,j] * p[group[i]+1]) 
        }
    }
    ## PRIORS 
    phi ~ dunif(0, 1)
    p[1] ~ dunif(0,1)
    p[2] ~ dunif(0,p[1])
    #one ~ dconstraint(p[1] >= p[2])   ## fixes lack of identifyablity
    pi ~ dunif(0, 1)
    
    for (i in 1:nind){
      for (t in 2:nyear){
        al[i,t-1] <- equals(z[i,t], 1)
      } #t
      for (t in 1:(nyear-1)){
        d[i,t] <- equals(z[i,t]-al[i,t],0)
      } #t   
      alive[i] <- sum(al[i,])
    } #i
    
    for (t in 1:(nyear-1)){
      N[t] <- sum(al[,t])        # Actual population size
      B[t] <- sum(d[,t])         # Number of entries
    } #t
}

## Form the list of data and constants
data <- list(nind = nind, 
                  nyear = nyear,
                  first = first,
                  y = y)

#data <- list(y = y, one = 1)

## Generate inits for the latent states:
z_init <- cbind(rep(NA, dim(z)[1]), z[,-1]) 

## Inits
inits <- function(){list(
    phi = 0.5,
    pi = 0.5,
    p = rep(0.5, 2),         
    group = rep(0, nind),
    z = z_init)}

## Run jags
het1 <- jags(model.file = code1, 
           data = data, 
           inits = inits,
           parameters.to.save = c("phi","pi","p","N"),
           n.iter  = 10000, 
           n.burnin = 5000, 
           n.chains = 5,
           n.thin=1)

het1$BUGSoutput


## Posterior inference
MCMCsummary(het1)

# mean         sd       2.5%       50%     97.5% Rhat n.eff
# p[1] 0.7548201 0.09765779 0.61419271 0.7312792 0.9739023 1.02    91
# p[2] 0.4263152 0.18168174 0.04754543 0.4627715 0.6757762 1.06    67
# phi  0.7140292 0.02843035 0.65961329 0.7134746 0.7707650 1.03   360
# pi   0.4400376 0.29269036 0.02104497 0.3778748 0.9558508 1.04    17

MCMCplot(het1)

#-- estimated

#-- truth
## p_class1: 0.7
## p_class2: 0.4
## prop_class1: 0.7   (note this is 1-pi)
## phi: 0.7


##--- 2bis. Capture-recapture model with constant parameters 
##---    and heterogeneity in detection using 2-class finite mixtures
##--- HMM formulation

## Model fitting
code2 <- function(){
  #DEFINE PARAMETERS
  #probabilities for each INITIAL STATES
  px0[1] <- pi # prob. of being in class 1
  px0[2] <- 1 - pi # prob. of being in class 2
  px0[3] <- 0 # prob. of being in initial state dead

  #OBSERVATION PROCESS: probabilities of observations (columns) at a given occasion given states (rows) at this occasion
  po[1,1] <- 1 - p[1]
  po[1,2] <- p[1]
  po[2,1] <- 1 - p[2]
  po[2,2] <- p[2]
  po[3,1] <- 1
  po[3,2] <- 0

  po.init[1,1] <- 0
  po.init[1,2] <- 1
  po.init[2,1] <- 0
  po.init[2,2] <- 1
  po.init[3,1] <- 1
  po.init[3,2] <- 0

  #STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
  px[1,1] <- phi
  px[1,2] <- 0
  px[1,3] <- 1 - phi
  px[2,1] <- 0
  px[2,2] <- phi
  px[2,3] <- 1 - phi
  px[3,1] <- 0
  px[3,2] <- 0
  px[3,3] <- 1
  ##
  for(i in 1:nind) { # for each ind
    z[i,first[i]] ~ dcat(px0[1:3])
    y[i,first[i]] ~ dcat(po.init[z[i,first[i]],1:2])
    for(j in (first[i]+1):nyear) { # loop over time
      ## STATE EQUATIONS ##
      ## draw states at j given states at j-1
      z[i,j] ~ dcat(px[z[i,j-1],1:3])
      ## OBSERVATION EQUATIONS ##
      ## draw observations at j given states at j
      y[i,j] ~ dcat(po[z[i,j],1:2])   
    }
  }
  ## PRIORS 
  phi ~ dunif(0, 1)
  p[1] ~ dunif(0,1)
  p[2] ~ dunif(0,p[1])
  # for(i in 1:2) {
  #   p[i] ~ dunif(0, 1)
  # }
  # one ~ dconstraint(p[1] >= p[2])   ## fixes lack of identifyablity
  pi ~ dunif(0, 1)
}

## Form the list of data and constants
bugs.data <- list(nind = nind, 
                  nyear = nyear,
                  first = first,
                  y = y +1)


## Generate inits
z_init <- cbind(rep(NA, dim(z)[1]), z[,-1]) 
mask <- sample(nind)
z_init[mask<floor(nind/2),] <- 1
inits <- list(
  phi = 0.5,
  pi = 0.5,
  p = rep(0.5, 2),
  z = z_init)

## Run jags
het2 <- jags(model.file = code2, 
             data = bugs.data, 
             inits = inits,
             parameters.to.save = c("phi","pi","p"),
             n.iter  = 10000, 
             n.burnin = 5000, 
             n.chains = 4,
             n.thin=1)

## Posterior inference
MCMCsummary(het2)

# mean         sd      2.5%       50%     97.5% Rhat n.eff
# p[1] 0.7138156 0.05063152 0.6223869 0.7112780 0.8134009 1.01   952
# p[2] 0.6326916 0.03799028 0.5560687 0.6344159 0.7027105 1.03   809
# phi  0.7003973 0.02380870 0.6538799 0.7006556 0.7467755 1.00  1195
# pi   0.2565821 0.04029415 0.1814802 0.2551673 0.3403638 1.01  1100

MCMCplot(het2)

#-- estimated

#-- truth
## p_class1: 0.7
## p_class2: 0.4
## prop_class1: 0.7   (note this is 1-pi)
## phi: 0.7

