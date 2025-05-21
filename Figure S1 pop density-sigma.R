# Manuscript Figure S!
# resident eq. population density as a function of its variance. 


# Parameter values
w       <- 1
theta1  <- 8
K1      <- 1

g1 <- 1/theta1

# Sequence of sigma1 values
sigma1.seq <- seq(0, 2.5, length.out = 100)  
V1.seq <- sigma1.seq^2

# Calculate n.EQ1 for each variance
n.EQ1.seq <- sapply(V1.seq, function(V1) {
  Omega1 <- sqrt(w^2 + 2*V1)
  (Omega1 / w) * (K1 - g1 * V1)
})

# Plot
plot(sigma1.seq, n.EQ1.seq, type = "l", lwd = 2,
     xlab = expression(sigma[1]),
     ylab = "Resident equilibrium density (n.EQ1)",
     main = "Resident equilibrium density vs variance")
grid()





