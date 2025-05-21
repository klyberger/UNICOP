# Manuscript Figure 6 
# coexistence outcomes

# check that these colors are correct

## -----------------------------------------------------------
## 1.  Parameter values (unchanged) ---------------------------
## -----------------------------------------------------------
w       <- 1
theta1  <- 8
theta2  <- 8
K1      <- 1
K2      <- 1
sigma1  <- 2.5
sigma2  <- 2.5

## -----------------------------------------------------------
## 2.  Helper functions -------------------
## -----------------------------------------------------------
g1 <- 1/theta1
g2 <- 1/theta2

frontL.ineq <- function(pars){
  with(pars, {
    V1     <- sigma1^2
    V2     <- sigma2^2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega1 <- sqrt(w^2 + 2*V1)
    g1     <- 1/theta1
    g2     <- 1/theta2
    g1*(V1 + (2*Omega^3)/Omega1)  # z.  (species 1 invades)
  })
}

backL.ineq <- function(pars){
  with(pars, {
    V1     <- sigma1^2
    V2     <- sigma2^2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega2 <- sqrt(w^2 + 2*V2)
    g2     <- 1/theta2
    g2*(V2 + (2*Omega^3)/Omega2)  # z.  (species 2 invades)
  })
}

fitness.z1 <- function(z1, pars){
  with(pars, {
    V1  <- sigma1^2
    V2  <- sigma2^2
    g1  <- 1/theta1
    g2  <- 1/theta2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega1 <- sqrt(w^2 + 2*V1)
    Omega2 <- sqrt(w^2 + 2*V2)
    n.EQ2  <- (Omega2/w)*(K2 - g2*V2)
    f1   <- K1 - g1*(z1^2 + V1)
    alpha11    <- w/Omega1^2
    alpha.inter<- (w/Omega)*exp(-(z1^2)/(2*Omega^2))
    f1 - alpha11*0 - alpha.inter*n.EQ2
  })
}

fitness.z2 <- function(z2, pars){
  with(pars, {
    V1  <- sigma1^2
    V2  <- sigma2^2
    g1  <- 1/theta1
    g2  <- 1/theta2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega1 <- sqrt(w^2 + 2*V1)
    Omega2 <- sqrt(w^2 + 2*V2)
    n.EQ1  <- (Omega1/w)*(K1 - g1*V1)
    f2   <- K2 - g2*(z2^2 + V2)
    alpha22    <- w/Omega2^2
    alpha.inter<- (w/Omega)*exp(-(z2^2)/(2*Omega^2))
    f2 - alpha22*0 - alpha.inter*n.EQ1
  })
}

zEQ1.function <- function(pars){
  with(pars, {
    V1  <- sigma1^2
    V2  <- sigma2^2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega2 <- sqrt(w^2 + 2*V2)
    g1  <- 1/theta1
    g2  <- 1/theta2
    -Omega*sqrt(-2*log( ((2*g1*Omega^3)/Omega2)/(K2 - g2*V2) ))
  })
}

zEQ2.function <- function(pars){
  with(pars, {
    V1  <- sigma1^2
    V2  <- sigma2^2
    Omega  <- sqrt(w^2 + V1 + V2)
    Omega1 <- sqrt(w^2 + 2*V1)
    g1  <- 1/theta1
    g2  <- 1/theta2
    -Omega*sqrt(-2*log( ((2*g2*Omega^3)/Omega1)/(K1 - g1*V1) ))
  })
}

## -----------------------------------------------------------
## 3.  Create sigma grids and classify scenarios -------------
## -----------------------------------------------------------
n.          <- 100
maxsig1     <- sqrt(K1*theta1)
maxsig2     <- sqrt(K2*theta2)
sigma1.seq  <- seq(0, maxsig1, length.out = n.)
sigma2.seq  <- seq(0, maxsig2, length.out = n.)
scenario.mat<- matrix(NA, nrow = n., ncol = n.)

which.caseI          <- c('22','21','20','12','11','10','02','01','00')
which.caseII         <- c('21','20','11','10','01','00')
which.caseII.reverse <- c('12','02','11','01','10','00')
which.caseIII        <- c('11','10','01','00')
caseIV               <- 20          # fitness = 0

for(ii in 1:n.){
  for(jj in 1:n.){
    pars <- list(w=w,K1=K1,K2=K2,
                 sigma1=sigma1.seq[ii], sigma2=sigma2.seq[jj],
                 theta1=theta1, theta2=theta2)

    f1 <- fitness.z1(0, pars)
    f2 <- fitness.z2(0, pars)

    if(frontL.ineq(pars) < K1){
      if(backL.ineq(pars) < K2){                         # ---- CASE I
        z1eq <- zEQ1.function(pars);  z2eq <- zEQ2.function(pars)
        if(!is.na(z1eq) & !is.na(z2eq)){
          who.1 <- ifelse(fitness.z1(0,   pars) > 0, 1, 0) +
                   ifelse(fitness.z1(z1eq, pars) > 0, 1, 0)
          who.2 <- ifelse(fitness.z2(0,   pars) > 0, 1, 0) +
                   ifelse(fitness.z2(z2eq, pars) > 0, 1, 0)
          scenario.mat[ii,jj] <- which(which.caseI == paste(who.1,who.2,sep=''))
        }
      } else {                                           # ---- CASE II (front<K1)
        z2eq <- zEQ2.function(pars)
        if(!is.na(z2eq)){
          who.1 <- ifelse(fitness.z1(0,   pars) > 0, 1, 0)
          who.2 <- ifelse(fitness.z2(0,   pars) > 0, 1, 0) +
                   ifelse(fitness.z2(z2eq, pars) > 0, 1, 0)
          scenario.mat[ii,jj] <- 9 + which(which.caseII.reverse ==
                                           paste(who.1,who.2,sep=''))
        }
      }
    } else {
      if(backL.ineq(pars) < K2){                         # ---- CASE II (back<K2)
        z1eq <- zEQ1.function(pars)
        if(!is.na(z1eq)){
          who.1 <- ifelse(fitness.z1(0,   pars) > 0, 1, 0) +
                   ifelse(fitness.z1(z1eq, pars) > 0, 1, 0)
          who.2 <- ifelse(fitness.z2(0,   pars) > 0, 1, 0)
          scenario.mat[ii,jj] <- 9 + which(which.caseII ==
                                           paste(who.1,who.2,sep=''))
        }
      } else {                                           # ---- CASE III
        who.1 <- ifelse(fitness.z1(0, pars) > 0, 1, 0)
        who.2 <- ifelse(fitness.z2(0, pars) > 0, 1, 0)
        scenario.mat[ii,jj] <- 15 + which(which.caseIII ==
                                          paste(who.1,who.2,sep=''))
      }
    }
  }
}
## -----------------------------------------------------------
## 4.  Outcome figure ------------------------
## -----------------------------------------------------------

outcome.cols <- c(
  "#2E8B57", "#2E8B57", "#7B3294", "#2E8B57", "#A0D995",
  "#C2A5CF", "#7B3294", "#C2A5CF", "#D01C8B",
  "#2E8B57", "#7B3294", "#A0D995", "#C2A5CF", "#7B3294",
  "#D01C8B", "#2E8B57", "#7B3294", "#7B3294", "#D01C8B", "white"
)

# unique legend entries
color_description_map <- c(
  "#2E8B57" = "persistence",
  "#A0D995" = "almost persistence",
  "#7B3294" = "competitive exclusion",
  "#C2A5CF" = "almost competitive exclusion",
  "#D01C8B" = "priority effect",
  "white"   = "fitness = 0"
)

all_cols   <- c(
  "#2E8B57",  # persistence
  "#A0D995",  # almost persistence
  "#7B3294",  # competitive exclusion
  "#C2A5CF",  # almost competitive exclusion
  "#D01C8B",  # priority effect
  "white"     # fitness = 0
)

present_cols <- intersect(all_cols, unique(outcome.cols[scenario.mat]))
legend_labs <- unname(color_description_map[present_cols])
breaks <- seq(0.5, 20.5, by = 1)

pdf("2spp_bifurcation_outcome_col_v3.pdf", width=6, height=5)

par(mfrow = c(1, 1), mar   = c(5, 4, 4, 10))
image(sigma1.seq, sigma2.seq, scenario.mat, col = outcome.cols,
      xlab = expression(sigma[1]), ylab = expression(sigma[2]),breaks=breaks)
legend("topright", inset = c(-.61,0),   legend = legend_labs,
  fill   = present_cols, bty = "n", xpd = TRUE, title = "Outcome")

dev.off()
