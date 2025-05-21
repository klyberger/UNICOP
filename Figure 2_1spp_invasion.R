# Manuscript figure 2

# Variance figure single invader
# 1 is invader
# 2 is resident

## Setting the parameters
w <- 1 ## competition width
theta1 <- 8 ## width of intrinsic growth function = 1/g
theta2 <- 8
K1 <- 1 ## vector of intrinsic growth potentials #1.7 for all scenarios
K2 <- 1  #2 for all scenarios
sigma1 <- 2.5 ## vector of species trait standard deviations
sigma2 <- 2.4

#check existence in 1 species
if(K1 < sigma1^2/theta1) print("error sp1")
if(K2 < sigma2^2/theta2) print("error sp2")


## ---------------------------------------------------------------------
## Testing the bimodality inequality varying sigma1 and sigma2 (sp 1 invades sp 2 resident)
## ---------------------------------------------------------------------

#max sigma^2 should be K*theta
maxsig1 = sqrt(K1*theta1)
maxsig2 = sqrt(K2*theta2)
  
n. <- 100
sigma1.seq <- seq(0, maxsig1, length.out = n.)
sigma2.seq <- seq(0, maxsig2, length.out = n.)

ineq1 <- matrix(NA, ncol = n., nrow = n.)

for(ii in 1:n.){
	for(jj in 1:n.){
	  Omega. <- sqrt(w^2 + sigma1.seq[ii]^2 + sigma2.seq[jj]^2)
		ineq1[ii, jj] <- ifelse((K2-sigma2.seq[jj]^2/theta2) > (2*Omega.^3)/(theta2 *sqrt(w^2 + 2*sigma2.seq[jj]^2)), 2, 1) #typo need sqrt for omega_2
	}
}

## ---------------------------------------------------------------------
## Also testing the sign of fitness varying sigma1 and sigma2 (sp 1 invades sp 2 resident)
## ---------------------------------------------------------------------

ineq2 <- matrix(NA, ncol = n., nrow = n.)

#A(A) is above peaks (and valleys) and B(B) is below peaks (and valleys)

for(ii in 1:n.){
  for(jj in 1:n.){
    Omega. <- sqrt(w^2 + sigma1.seq[ii]^2 + sigma2.seq[jj]^2)
    ineq2[ii, jj] <- ifelse(
      (K2-sigma2.seq[jj]^2/theta2) > (2*Omega.^3)/(theta2 *sqrt(w^2 + 2*sigma2.seq[jj]^2)),#bimodal #typo need sqrt for omega_2
      ifelse(
        (K1-sigma1.seq[ii]^2/theta1)/(K2-sigma2.seq[jj]^2/theta2) > sqrt((w^2 + 2*sigma2.seq[jj]^2)/Omega.^2),#invades at all eq. #typo Omega^2
        "bimodal-persistence",
        ifelse(
          (K1-sigma1.seq[ii]^2/theta1) > 2*Omega.^2/theta2 * (1-log((2*Omega.^3)/(theta2 *sqrt(w^2 + 2*sigma2.seq[jj]^2)*(K2-sigma2.seq[jj]^2/theta2)))), # invades only at peaks #typo need sqrt for omega_2
          "bimodal-almost persistence",
          "bimodal-no invasion" #no invasion
        )
      ),
      ifelse(
        (K1-sigma1.seq[ii]^2/theta1)/(K2-sigma2.seq[jj]^2/theta2) > sqrt((w^2 + 2*sigma2.seq[jj]^2)/Omega.^2), #typo Omega^2
        "unimodal-persistence", #unimodal, invades at 0
        "unimodal-no invion" #unimodal, doesn't invade at 0
      )
    )
  }
}

# Define the mapping of categories to numbers
category_map <- c("bimodal-persistence" = 1, "bimodal-almost persistence" = 2, "bimodal-no invasion" = 3, "unimodal-persistence" = 4, "unimodal-no invion" = 5)

# Create a numeric version of ineq2
ineq2_numeric <- matrix(category_map[as.character(ineq2)], nrow = nrow(ineq2))

# Define the color palette
my.colors <- c("#ff6e00", "#ffaa00", "#ffc800", "#6464ff", "#c8c8ff")


### Figure 2 ####
pdf(file="1spp_bifurcation_v3.pdf", width=6, height=4.7)


# Create the image
par(mar = c(5, 4.3, 4, 12))
#par(mfrow = c(1,1))
image(sigma1.seq, sigma2.seq, ineq2_numeric, col = my.colors, xlab=expression(V[1]^2), ylab=expression(V[2]^2))

my.colors <- setNames(
  c("#ff6e00", "#ffaa00", "#ffc800", "#6464ff", "#c8c8ff"),
  names(category_map)
)
present <- intersect(
  names(category_map),
  unique(as.character(ineq2))      
)

# Add a legend for those
legend(x = "topright",inset  = c(-0.9, 0),legend = present,fill   = my.colors[present],title  = "Case",bty    = "n",bg     = NA, xpd    = TRUE
)



## ---------------------------------------------------------------------
## Example fitness functions to overlay
## ---------------------------------------------------------------------
# set up cases A, B, D, lowerright_D, E (no case C found yet)
sigma1_vec <- c(1, 0.1, 1.5, 1.5, 2.3)
sigma2_vec <- c(0.1, 1, 0.5, 2, 0.1)

z1. <- seq(-3, 3, by = 0.1)

# shared parameters list
pars_base <- list(w = w, K1 = K1, K2 = K2, theta1 = theta1, theta2 = theta2)


# Define the function
r. <- function(z1, pars){
    # read in the parameters
    sigma1 <- pars$sigma1
    sigma2 <- pars$sigma2
    w <- pars$w
    K1 <- pars$K1
    K2 <- pars$K2
    theta1 <- pars$theta1
    theta2 <- pars$theta2

    # Calculations
    n.EQ1 <- 0 # We want this one to be zero
    n.EQ2 <- sqrt(1 + 2*sigma2^2/w^2)*(K2 - (1/theta2)*sigma2^2)
    Omega. <- sqrt(w^2 + sigma1^2 + sigma2^2) 

    dz <- z1-0
    b1 <- K1 - (z1^2 + sigma1^2)/theta1
    alpha1.0 <- 1*w/sqrt(sigma1^2+sigma1^2+w^2)
    alpha.inter <- exp(-dz^2/(2*(sigma1^2+sigma2^2+w^2)))*w/sqrt(sigma1^2+sigma2^2+w^2) 

    out. <- b1 - alpha1.0*n.EQ1 - alpha.inter*n.EQ2
    return(out.)
}

# --- Case A inset (shifted left by 0.05) ---
par(fig = c(0.2, 0.37, 0.18, 0.38), new = TRUE, mar = c(1, 1, 1, 1))
pars <- modifyList(pars_base, list(sigma1 = sigma1_vec[1], sigma2 = sigma2_vec[1]))
plot(z1., r.(z1., pars), type='l', xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(-0.4,0.2))
abline(h=0, lty=2)

# --- Case B inset ---
par(fig = c(0.12, 0.29, 0.32, 0.52), new = TRUE, mar = c(1, 1, 1, 1))
pars <- modifyList(pars_base, list(sigma1 = sigma1_vec[2], sigma2 = sigma2_vec[2]))
plot(z1., r.(z1., pars), type='l', xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(-0.4,0.2))
abline(h=0, lty=2)

# --- Case D inset ---
par(fig = c(0.28, 0.45, 0.60, 0.80), new = TRUE, mar = c(1, 1, 1, 1))
pars <- modifyList(pars_base, list(sigma1 = sigma1_vec[4], sigma2 = sigma2_vec[4]))
plot(z1., r.(z1., pars), type='l', xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(-0.4,0.2))
abline(h=0, lty=2)

# --- Case E inset ---
par(fig = c(0.45, 0.62, 0.40, 0.60), new = TRUE, mar = c(1, 1, 1, 1))
pars <- modifyList(pars_base, list(sigma1 = sigma1_vec[5], sigma2 = sigma2_vec[5]))
plot(z1., r.(z1., pars), type='l', xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(-0.4,0.2))
abline(h=0, lty=2)


dev.off()
