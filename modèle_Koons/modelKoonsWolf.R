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

#--------------------------------------------------------
# 2. Process model of state-space likelihood
#--------------------------------------------------------
for (t in 1:years) {
  N[t+1] <- N[t] + eps[t]
  eps[t] ~ dnorm(0,tauPro)

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
Yse=rep(0.3,27)

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

################################################################################
#Attention code ci-dessous en format Rmd
################################################################################

## Modèle Marche aléatoires


$$N[t+1]=N[t]+eps$$
  $$eps \sim \text{Normale}(0,\sigma_\text{proc})$$
  $$\sigma_\text{proc} \sim \text{Uniforme}(0,5)$$
  $$y[t] \sim \text{Normale}(N[t],\text{ySe}[t]^{-2})$$
  ```{r}
modelnull = function(){
  # Priors 
  sigmaProc ~ dunif(0,5)
  tauProc = pow(sigmaProc,-2)
  
  N1 ~ dnorm(17,1/0.28^2)
  N[1] = max(0,N1)
  
  # Process model
  for (t in 1:years){
    N[t+1] = N[t] + eps[t]
    eps[t] ~ dnorm(0,tauProc)
  }
  
  # Observation model
  for (t in 1:years){
    y[t] ~ dnorm(N[t],pow(yse[t],-2))
  }
}

```



```{r,echo=FALSE}
perspredictwolf=matrix(NA,27,1)
perspredictwolf[1]=dat$N[1]

quantile=matrix(NA,27,2)
quantile[1,1] = quantile(dat$N[1], probs = 2.5/100)
quantile[1,2] = quantile(dat$N[1], probs = 97.5/100)

for (i in 2:length(dat$N)){
  y = dat$N[1:i]
  yse = ObsSE[1:i]
  years = i
  # Initialisation des données
  bugs.data = list(
    years = years,
    y = y,
    yse = yse) 
  
  
  # Paramètres JAGS
  bugs.monitor = c("sigmaProc","N", "tauProc")
  bugs.chains = 3
  bugs.inits = function(){
    list(
    ) 
  }
  
  # Lancement du programme
  wolf_modelnull = jags(data = bugs.data,
                        inits = bugs.inits, 
                        parameters.to.save = bugs.monitor,
                        model.file = modelnull,
                        n.chains = bugs.chains, 
                        n.thin=10,
                        n.iter=20000, 
                        n.burnin=5000)
  
  perspredictwolf[i] = median(wolf_modelnull$BUGSoutput$sims.list$N[,i+1])
  quantile[i,1] = quantile(wolf_modelnull$BUGSoutput$sims.list$N[,i+1], probs = 2.5/100)
  quantile[i,2] = quantile(wolf_modelnull$BUGSoutput$sims.list$N[,i+1], probs = 97.5/100)
  
}
```

```{r}
ggplot()+
  geom_line(aes(x = 1995 + 1:unique(nyears), y = perspredictwolf), colour = "red", lty = "dashed")+
  geom_ribbon(aes(x = 1995 + 1:unique(nyears), ymin = quantile[,1], ymax = quantile[,2]), fill = "red", alpha = 0.3)+
  geom_point(data = bugs.data %>% as_tibble, aes(x = 1995 + 1:unique(nyears), y = y)) + 
  #coord_cartesian(xlim=c(1996,2023),ylim=c(0,1250))+
  theme_bw()+
  labs(title = "Estimated population size",
       x = "Years",
       y = "Number of wolves")
```

### Prédiction
Initialisation des différents taux de prélèvement :
  ```{r}
dH = c(0, 0.10, 0.20, 0.30)
```

```{r}
modelnull = function(){
  # Priors 
  sigmaProc ~ dunif(0,5)
  tauProc = pow(sigmaProc,-2)
  
  N1 ~ dnorm(17,1/0.28^2)
  N[1] = max(0,N1)
  
  # Process model
  for (t in 1:(years)){
    N[t+1] = N[t] + eps[t]
    eps[t] ~ dnorm(0,tauProc)
  }
  
  # Observation model
  for (t in 1:years){
    y[t] ~ dnorm(N[t],pow(yse[t],-2))
  }
  
  # Projection model
  for (t in (years+1):(years+2)){
    N[t+1] = (1-h)*N[t] + eps[t]
    eps[t] ~ dnorm(0,tauProc)
  }
  
}

```


```{r}
perspredictwolf=matrix(NA,29,1)
perspredictwolf[1]=dat$N[1]

quantile=matrix(NA,29,2)
quantile[1,1] = quantile(dat$N[1], probs = 2.5/100)
quantile[1,2] = quantile(dat$N[1], probs = 97.5/100)

ObsSE[28]=0.28
ObsSE[29]=0.3

for (i in 2:(length(dat$N)+2)){
  y = dat$N[1:i]
  yse = ObsSE[1:i]
  years = i
  h=0
  if (i>length(dat$N)){h=dH[1]}
  
  
  # Initialisation des données
  bugs.data = list(
    years = years,
    y = y,
    yse = yse,
    h = h) 
  
  
  # Paramètres JAGS
  bugs.monitor = c("sigmaProc","N", "tauProc")
  bugs.chains = 3
  bugs.inits = function(){
    list(
    ) 
  }
  
  # Lancement du programme
  wolf_modelnull = jags(data = bugs.data,
                        inits = bugs.inits, 
                        parameters.to.save = bugs.monitor,
                        model.file = modelnull,
                        n.chains = bugs.chains, 
                        n.thin=10,
                        n.iter=20000, 
                        n.burnin=5000)
  
  perspredictwolf[i] = median(wolf_modelnull$BUGSoutput$sims.list$N[,i+1])
  quantile[i,1] = quantile(wolf_modelnull$BUGSoutput$sims.list$N[,i+1], probs = 2.5/100)
  quantile[i,2] = quantile(wolf_modelnull$BUGSoutput$sims.list$N[,i+1], probs = 97.5/100)
  
}
```

```{r,warning=FALSE}
ggplot()+
  geom_line(aes(x = 1995 + 1:(unique(nyears)+2), y = perspredictwolf), colour = "red", lty = "dashed")+
  geom_ribbon(aes(x = 1995 + 1:(unique(nyears)+2), ymin = quantile[,1], ymax = quantile[,2]), fill = "red", alpha = 0.3)+
  geom_point(data = bugs.data %>% as_tibble, aes(x = 1995 + 1:(unique(nyears)+2), y = y)) + 
  coord_cartesian(xlim=c(1996,2023),ylim=c(0,1250))+
  theme_bw()+
  labs(title = "Projected and estimated population size",
       subtitle = "Harvest rate : 0%",
       x = "Years",
       y = "Number of wolves")
```


