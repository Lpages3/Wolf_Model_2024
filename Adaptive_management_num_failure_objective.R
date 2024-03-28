library(R2jags)
library(tidyverse)

meanfail1 = c()
meanfail2 = c()

T1 = Sys.time()

for (n in 1:15) {
  nyears = 25
  N1 = 30
  
  ssm_sim4 = data.frame(
    Year = 1:nyears,
    y = numeric(nyears),
    N = numeric(nyears),
    ybis = numeric(nyears),
    Nbis = numeric(nyears)
  )
  
  ssm_sim4$N[1] = N1
  ssm_sim4$Nbis[1] = N1
  # Paramètres initiaux
  pas = 1
  H = 0
  sigma = 0.15
  K = 800
  alpha = 0.5
  ite = 0
  tempH = c()
  NAMharvest = 0.15
  
  # Lancement du modèle
  for (nyears in seq(2, nyears, pas)) {
    # Boucle sur le nombre de tranches d'années parcourues
    ite = ite + 1
    tempH[ite] = H # Enregistre les taux de prélevement pour chaque année
    
    # Simulation des effectifs
    if (nyears <= 5) {
      # Initialisation des effectifs, sans prélevement sur les 5 premières années
      for (t in 1:(nyears - 1)) {
        u = ssm_sim4$N[t]
        
        Er = exp(alpha * (1 - u / K)) * u
        
        ssm_sim4$N[t + 1] = rpois(1, Er)
        ssm_sim4$Nbis[t + 1] = rpois(1, Er)
      }
    }
    
    if (nyears > 5) {
      # Suite de la simulation des effectifs entre les années 6 et 25
      for (t in (nyears - pas):(nyears - 1)) {
        u = ssm_sim4$N[t] * (1 - H)             # Taux de prélevement adaptatif
        v = ssm_sim4$Nbis[t] * (1 - NAMharvest) # Taux de prélevement constant
        
        Er = exp(alpha * (1 - c(u, v) / K)) * c(u, v)
        
        ssm_sim4$N[t + 1] = rpois(1, Er[1])
        ssm_sim4$Nbis[t + 1] = rpois(1, Er[2])
      }
    }
    # Simulation des effecitfs observés
    for (t in 1:nyears) {
      ssm_sim4$y[t] = rpois(1, ssm_sim4$N[t])
      ssm_sim4$ybis[t] = rpois(1, ssm_sim4$Nbis[t])
    }
    
    # Début de l'estimation par approche bayésienne
    # Initialisation des données
    bugs.data = list(nyears = nyears,
                     y = c(ssm_sim4$y[1:nyears], rep(NA, 5)),
                     dH = H)
    
    # Paramètres JAGS
    bugs.monitor = c("alpha", "sigmaProc", "tauProc", "K", "N")
    bugs.chains = 3
    init1 = list(alpha = .5, sigmaProc = .25)
    init2 = list(alpha = .1, sigmaProc = .05)
    init3 = list(alpha = 1, sigmaProc = .45)
    bugs.inits = list(init1, init2, init3)
    
    # Lancement du modèle
    
    wolf_modellogist = jags(
      data = bugs.data,
      inits = bugs.inits,
      parameters.to.save = bugs.monitor,
      model.file = modellogist,
      n.chains = bugs.chains,
      n.thin = 10,
      n.iter = 20000,
      n.burnin = 5000
    )
    output1 = wolf_modellogist$BUGSoutput$sims.matrix
    
    # Calcul du taux de reproduction estimé
    Nest = wolf_modellogist$BUGSoutput$median$N
    l = length(Nest)
    lamb = c()
    for (t in 1:(l - 5)) {
      lamb[t] = Nest[t + 1] / Nest[t]
    }
    lambda = mean(lamb)
    print(lambda)
    
    # Conditions de modification du taux de prélevement
    if (lambda < 1.2) {
      H = 0
    }
    if (lambda >= 1.2 & lambda < 1.3) {
      H = 0.1
    }
    if (lambda >= 1.3 & lambda < 1.4) {
      H = 0.2
    }
    if (lambda > 1.4) {
      H = 0.3
    }
    
    print(H)
  }
  
  
  
  
  # Initialisation des données
  bugs.data = list(nyears = nyears,
                   y = c(ssm_sim4$ybis[1:nyears], rep(NA, 5)),
                   dH = NAMharvest)
  
  wolf_modellogist = jags(
    data = bugs.data,
    inits = bugs.inits,
    parameters.to.save = bugs.monitor,
    model.file = modellogist,
    n.chains = bugs.chains,
    n.thin = 10,
    n.iter = 20000,
    n.burnin = 5000
  )
  
  output2 = wolf_modellogist$BUGSoutput$sims.matrix %>%
    as_tibble() %>%
    pivot_longer(cols = everything(),
                 values_to = "value",
                 names_to = "parameter2") %>%
    filter(str_detect(parameter2, "N")) %>%
    group_by(parameter2) %>%
    summarize(
      medianN2 = median(value),
      lq2 = quantile(value, probs = 2.5 / 100),
      hq2 = quantile(value, probs = 97.5 / 100)
    ) %>%
    mutate(years = parse_number(parameter2)) %>%
    arrange(years)
  
  output1 = output1 %>%
    as_tibble() %>%
    pivot_longer(cols = everything(),
                 values_to = "value",
                 names_to = "parameter1") %>%
    filter(str_detect(parameter1, "N")) %>%
    group_by(parameter1) %>%
    summarize(
      medianN1 = median(value),
      lq1 = quantile(value, probs = 2.5 / 100),
      hq1 = quantile(value, probs = 97.5 / 100)
    ) %>%
    mutate(years = parse_number(parameter1)) %>%
    arrange(years)
  
  
  output = output1 %>% left_join(output2)
  
  
  # Ranger les taux de croissances de chaque méthode
  lambda1 = c()
  lambda2 = c()
  
  for (t in 1:(nyears - 1)) {
    lambda1[t] = output1$medianN1[t + 1] / output1$medianN1[t]
    lambda2[t] = output2$medianN2[t + 1] / output2$medianN2[t]
  }
  
  # Compter le nombre de fois que le taux de croissance n'est pas dans l'objectif [1,1.2]
  fail1 = 0
  fail2 = 0
  for (t in 1:length(lambda1)) {
    if (between(lambda1[t], 1, 1.2) == FALSE) {
      fail1 = fail1 + 1
    }
    if (between(lambda2[t], 1, 1.2) == FALSE) {
      fail2 = fail2 + 1
    }
  }
  
  meanfail1[n] = fail1
  meanfail2[n] = fail2
}

mean(meanfail1)/nyears
mean(meanfail2)/nyears

temps = (Sys.time() - T1)/60

print (c("Le programme a duré :", Sys.time() - T1, "minutes."))