# Manuscript figure 6B
# Temporal dynamics 

## ---------------------------------------------------------------------
## Temporal Dynamics
## ---------------------------------------------------------------------
require(deSolve)
require(tidyverse)
require(cowplot)
require(dplyr)

# case 1AA (Persistence), 3EE (Competitive Exclusion), 3DE (Priority Effects)

organize_results <- function(sol, pars) {
  dat <- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  names(dat)[1] <- "time" ## name the first column "time"
  names(dat)[2:3] <- paste0("n_", 1:2) ## name abundance columns (n_k)
  names(dat)[4:5] <- paste0("m_", 1:2) ## name trait mean columns (m_k)
  dat <- dat %>%
    gather("variable", "v", 2:ncol(dat)) %>% ## long format
    separate(variable, c("type", "species"), sep="_") %>%
    spread(type, v) %>% ## separate columns for densities n and trait means m
    select(time, species, n, m) %>% ## rearrange columns
    mutate(species=as.integer(species), sigma1=pars$sigma1, w=pars$w,
           theta=pars$theta) %>% ## add parameters
    mutate(species=as.character(species))
  return(dat)
}

## cuttoff function from pastore:
## - array of values with smoothed step function applied to them
cutoff <- function(n) {
  return(ifelse(n<1, (1*(n>0))*(n*n*n*(10+n*(-15+6*n))), 1))
}

eqs <- function(time, state, pars) {
  n1 <- state[1] ## species densities
  n2 <- state[2]
  z1 <- state[3] ## species trait means
  z2 <- state[4]
  
  b1 <- K1 - (z1^2 + sigma1^2)/theta1
  b2 <- K2 - (z2^2 + sigma2^2)/theta2 #eqn S32 from Pastore suppl
  
  g1 <- -2*z1*sigma1^2/theta1
  g2 <- -2*z2*sigma2^2/theta2 #eqn S33 from Pastore suppl
  dz <- z1-z2 ## difference matrix of trait means
  
  alpha1.0 <- 1*w/sqrt(sigma1^2+sigma1^2+w^2)
  alpha2.0 <- 1*w/sqrt(sigma2^2+sigma2^2+w^2) #eqn S34 from Pastore suppl
  alpha.inter <- exp(-dz^2/(2*(sigma1^2+sigma2^2+w^2)))*w/sqrt(sigma1^2+sigma2^2+w^2)
  
  beta1.0 <- 0
  beta2.0 <- 0
  beta1 <- alpha.inter*2*sigma1^2*(-dz)/(2*(sigma1^2+sigma2^2+w^2)) #eqn S35 pastore suppl **check these**
  beta2 <- alpha.inter*2*sigma2^2*(dz)/(2*(sigma1^2+sigma2^2+w^2))
  ## The dndt equations are multiplied by a cutoff function, to make
  ## behavior smooth as abundances approach 0. The cutoff function is
  ## a smoothed-out step function whose derivative exists everywhere.
  dn1dt <- n1*(b1 - alpha1.0*n1 - alpha.inter*n2)*cutoff(n1/(1e-7)) #eqn S36 pastore suppl
  dn2dt <- n2*(b2 - alpha2.0*n2 - alpha.inter*n1)*cutoff(n2/(1e-7))
  dz1dt <- (g1 - beta1.0*n1 - beta1*n2) #removed h2 from front
  dz2dt <- (g2 - beta2.0*n2 - beta2*n1) #eqn S31 pastore suppl
  ## return equations by first flattening them back into a single vector
  return(list(c(dn1dt, dn2dt, dz1dt, dz2dt)))
}

## Listing the parameters
## case 1AA (Persistence)
w <- 1 ## competition width #3
theta1 <- 8 ## width of intrinsic growth function #8
theta2 <- 8
K1 <- 1 
K2 <- 1 ## vector of intrinsic growth potentials
sigma1 <- 1 ## vector of species trait standard deviations #2.5
sigma2 <- 0.5

n.EQ1 <- sqrt(1 + 2*sigma1^2/w^2)*(K2 - (1/theta1)*sigma1^2)
n.EQ2 <- sqrt(1 + 2*sigma2^2/w^2)*(K2 - (1/theta2)*sigma2^2)

tmax <- 250 ## time to integrate equations for
stepout <- tmax/250 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

ninit <- c(0.01,n.EQ2) ## initial densities 
muinit <- c(2,0) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

pars <- list(w = w, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

sol <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

ninit <- c(n.EQ1, 0.01) ## initial densities 
muinit <- c(0,2) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

sol2 <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

## create plot of species densities through time
densplot <- ggplot() +
  labs(title = "    Persistence") +
  geom_line(data=sol, aes(x=time, y=n, colour=as.factor(species))) +
  geom_line(data=sol2, aes(x=time, y=n, colour=as.factor(species)), linetype="dashed") +
  scale_y_continuous(name="pop. density", limits=c(0, NA)) +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=3, fontface="bold") +
  theme(legend.position="none")

## create plot of species traits through time
traitplot <- ggplot() +
  # geom_ribbon(data = sol, aes(x=time, ymin=m-sigma1, ymax=m+sigma1,
  #                 fill=as.factor(species)), alpha=0.15) +
  geom_line(data = sol, aes(x=time, y=m, colour=as.factor(species))) +
  geom_line(data = sol2, aes(x=time, y=m, colour=as.factor(species)),linetype="dashed") +
  ylab("trait value") +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  scale_fill_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=0.70, fontface="bold") +
  theme(legend.position="none")

plot_grid(densplot, traitplot, ncol = 1, nrow = 2) #solid spp1 (blue) invades, dashed spp2 (yellow) invades

## Case 3DD (Priority Effect)
sigma1 <- 2 ## vector of species trait standard deviations 
sigma2 <- 2

n.EQ1 <- sqrt(1 + 2*sigma1^2/w^2)*(K2 - (1/theta1)*sigma1^2)
n.EQ2 <- sqrt(1 + 2*sigma2^2/w^2)*(K2 - (1/theta2)*sigma2^2)

tmax <- 250 ## time to integrate equations for
stepout <- tmax/250 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

ninit <- c(0.01,n.EQ2) ## initial densities 
muinit <- c(2,0) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

pars <- list(w = w, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

sol <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

ninit <- c(n.EQ1, 0.01) ## initial densities 
muinit <- c(0,2) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

sol2 <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

## create plot of species densities through time
densplotPE <- ggplot() +
  labs(title = "    Priority Effect") +
  geom_line(data=sol, aes(x=time, y=n, colour=as.factor(species))) +
  geom_line(data=sol2, aes(x=time, y=n, colour=as.factor(species)), linetype="dashed") +
  scale_y_continuous(name="pop. density", limits=c(0, NA)) +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=3, fontface="bold") +
  theme(legend.position="none")

## create plot of species traits through time
traitplotPE <- ggplot() +
  # geom_ribbon(data = sol, aes(x=time, ymin=m-sigma1, ymax=m+sigma1,
  #                 fill=as.factor(species)), alpha=0.15) +
  geom_line(data = sol, aes(x=time, y=m, colour=as.factor(species))) +
  geom_line(data = sol2, aes(x=time, y=m, colour=as.factor(species)),linetype="dashed") +
  ylab("trait value") +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  scale_fill_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=0.70, fontface="bold") +
  theme(legend.position="none")

plot_grid(densplotPE, traitplotPE, ncol = 1, nrow = 2) #solid spp1 (blue) invades, dashed spp2 (yellow) invades *but this looks like PE?

## Case 3EE (Permanence)
sigma1 <- 1.8 ## vector of species trait standard deviations 
sigma2 <- 0.5

## Case 3DE (Competitive Exclusion)
sigma1 <- 0.5 ## vector of species trait standard deviations 
sigma2 <- 2.5

n.EQ1 <- sqrt(1 + 2*sigma1^2/w^2)*(K2 - (1/theta1)*sigma1^2)
n.EQ2 <- sqrt(1 + 2*sigma2^2/w^2)*(K2 - (1/theta2)*sigma2^2)

tmax <- 250 ## time to integrate equations for
stepout <- tmax/250 ## time step size for output
time <- seq(0, tmax, by=stepout) ## sampling points in time

ninit <- c(0.01,n.EQ2) ## initial densities 
muinit <- c(2,0) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

pars <- list(w = w, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

sol <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

ninit <- c(n.EQ1, 0.01) ## initial densities 
muinit <- c(0,2) ## initial trait means 
ic <- c(ninit, muinit) ## initial conditions coerced into a vector

sol2 <- ode(func=eqs, y=ic, parms=pars, times=time) %>% ## solve ODEs
  organize_results(pars) ## put results in tidy table

## create plot of species densities through time
densplotCE <- ggplot() +
  labs(title = "    Competitive Exclusion") +
  geom_line(data=sol, aes(x=time, y=n, colour=as.factor(species))) +
  geom_line(data=sol2, aes(x=time, y=n, colour=as.factor(species)), linetype="dashed") +
  scale_y_continuous(name="pop. density", limits=c(0, NA)) +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=3, fontface="bold") +
  theme(legend.position="none")

## create plot of species traits through time
traitplotCE <- ggplot() +
  # geom_ribbon(data = sol, aes(x=time, ymin=m-sigma1, ymax=m+sigma1,
  #                 fill=as.factor(species)), alpha=0.15) +
  geom_line(data = sol, aes(x=time, y=m, colour=as.factor(species))) +
  geom_line(data = sol2, aes(x=time, y=m, colour=as.factor(species)),linetype="dashed") +
  ylab("trait value") +
  scale_colour_manual(values=c("#0072B2", "#E69F00")) +
  scale_fill_manual(values=c("#0072B2", "#E69F00")) +
  annotate(geom="text", label="", x=0, y=0.70, fontface="bold") +
  theme(legend.position="none")

plot_grid(densplotCE, traitplotCE, ncol = 1, nrow = 2) #solid spp1 (blue) invades, dashed spp2 (yellow) invades

combinedPlot1 <- plot_grid(densplot, traitplot, labels = "1AA", ncol = 1)
combinedPlot2 <- plot_grid(densplotPE, traitplotPE, labels = "3DD", ncol = 1)
combinedPlot3 <- plot_grid(densplotCE, traitplotCE, labels = "3DE", ncol = 1)

# Now, combine all the plots side by side
finalCombinedPlot <- plot_grid(combinedPlot1, combinedPlot2, combinedPlot3, ncol = 3, align = 'v')

# Export as PDF
pdf("combined_plots.pdf", width = 10, height = 7)
print(finalCombinedPlot)
dev.off()

