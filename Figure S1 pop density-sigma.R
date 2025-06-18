# Manuscript Figure S!
# resident eq. population density as a function of its variance. 


# Parameter values
w       <- 1
theta1  <- 8
K1      <- 1
g1 <- 1/theta1

# Sequence of sigma1 values
V1.seq <- seq(0, 6.25, length.out = 100)

# Calculate n.EQ1 for each variance
n.EQ1.seq <- sapply(V1.seq, function(V1) {
  Omega1 <- sqrt(w^2 + 2*V1)
  (Omega1 / w) * (K1 - g1 * V1)
})

pdf("FigureS1_density-variance.pdf", width = 6, height = 4)

# Plot
plot(V1.seq, n.EQ1.seq, type = "l", lwd = 2,
     xlab = expression(V[1]),
     ylab = "Resident equilibrium density",
     main = "")
grid()

dev.off()


# Parameter ranges
w.values  <- c(1, 5)
K1.values <- c(0.5, 1, 5)
g1.values <- c(0.125, 1)

# Sequence of variance values (evenly spaced)
V1.seq <- seq(0, 6.25, length.out = 100)

# Set up plotting grid: 2 (w) × 3 (K1) × 2 (g1) = 12 plots
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))  # 4 rows × 3 columns

for (w in w.values) {
  for (K1 in K1.values) {
    for (g1 in g1.values) {
      
      # Compute equilibrium densities
      n.EQ1.seq <- sapply(V1.seq, function(V1) {
        Omega1 <- sqrt(w^2 + 2 * V1)
        (Omega1 / w) * (K1 - g1 * V1)
      })
      
      # Plot
      plot(V1.seq, n.EQ1.seq, type = "l", lwd = 2, col = "darkblue",
           xlab = expression(V[1]),
           ylab = "n.EQ1", 
           main = bquote(w == .(w) ~ ", " ~ K[1] == .(K1) ~ ", " ~ g[1] == .(g1)))
      abline(h = 0, lty = 3, col = "gray")  # zero line
      grid()
    }
  }
}



