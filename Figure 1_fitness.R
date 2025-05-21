#Manuscript Figure 1 

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


## Setting the parameters
w <- 1 ## competition width
theta1 <- 8 ## width of intrinsic growth function = 1/g
theta2 <- 8
K1 <- 2 #** increase K
K2 <- 1 ## vector of intrinsic growth potentials
sigma1 <- 4 ## vector of species trait standard deviations
sigma2 <- 4

z1. <- seq(-10, 10, by = 0.1)

#fitness function across z

#for high/low K
sigma1 < sqrt(5*theta1) #check species can persist alone
sigma2 < sqrt(K2*theta2)

pars.KH <- list(w = w, K1 = 5, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

pars.KL <- list(w = w, K1 = 0.5, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

#for high/low g = 1/theta
sigma1 < sqrt(K1*10) #check species can persist alone

pars.gH <- list(w = w, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = 1, theta2 = 1)

pars.gL <- list(w = w, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = 10, theta2 = 10)

#for high/low V = sigma^2
pars.VH <- list(w = w, K1 = K1, K2 = K2, sigma1 = 6, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

pars.VL <- list(w = w, K1 = K1, K2 = K2, sigma1 = 0.5, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)



#alpha function across z

# Define the function
a. <- function(z1, pars){
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
    return(alpha.inter)
}

#for high/low w
pars.wH <- list(w = 5, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

pars.wL <- list(w = 1, K1 = K1, K2 = K2, sigma1 = sigma1, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

#for high/low V = sigma^2
6 < sqrt(K1*theta1) #check species can persist alone

pars.V1H <- list(w = w, K1 = K1, K2 = K2, sigma1 = 6, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

pars.V1L <- list(w = w, K1 = K1, K2 = K2, sigma1 = 0.5, sigma2 = sigma2, theta1 = theta1, theta2 = theta2)

### create figure ####
par(mfrow=c(1, 2), mar=c(4,4,2,0), oma=c(1,1,1,1))
pdf(file="Fitness_2_panel.pdf", width=5.5, height=3.5)
layout(matrix(c(1, 1, 1, 2, 2), ncol=5, byrow=TRUE))
plot(z1.,r.(z1 = z1., pars.KH), type = 'l', xlab="invader trait", ylab="fitness", ylim=c(-3,5))
lines(z1.,r.(z1 = z1., pars.KL), type = 'l',lty ="dashed")
legend(x=4, y=5.5, legend=c(expression(paste(K[1], " = 5")), expression(paste(K[1], " = 0.55")),"g = 1", "g = 0.1",expression(paste(V[1], " = 36")), expression(paste(V[1], " = 0.25"))), col=c("black", "black","blue", "blue", "red", "red"), lty=c(1,2), bty="n",seg.len=1, y.intersp=1, x.intersp=0.75)
lines(z1.,r.(z1 = z1., pars.gH), type = 'l', col="blue")
lines(z1.,r.(z1 = z1., pars.gL), type = 'l', lty ="dashed", col="blue")
lines(z1.,r.(z1 = z1., pars.VH), type = 'l', col = "red")
lines(z1.,r.(z1 = z1., pars.VL), type = 'l', lty ="dashed", col="red")

#alpha panel
plot(z1.,a.(z1 = z1., pars.wH), type = 'l', xlab=expression(paste(Delta, "z")), ylab="alpha", ylim=c(0,1), col="dark green")
lines(z1.,a.(z1 = z1., pars.wL), type = 'l',col="dark green", lty="dashed")
legend(x=-1.25, y=1.05, legend=c(expression(paste(omega," = 5")), expression(paste(omega," = 1")), expression(paste(V[1], " = 36")), expression(paste(V[1], " = 0.25"))), col=c("dark green", "dark green","red", "red"), lty=c(1,2), bty="n",seg.len=1, y.intersp=1, x.intersp=0.75)
lines(z1.,a.(z1 = z1., pars.V1H), type = 'l', ylab="", ylim=c(0,1),col="red")
lines(z1.,a.(z1 = z1., pars.V1L), type = 'l',col="red", lty = "dashed")
dev.off()

