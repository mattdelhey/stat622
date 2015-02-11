set.seed(1)

data <- c(2.3656491, 2.4952035, 1.0837817, 0.7586751, 0.8780483, 1.2765341,
          1.4598699, 0.1801679, -1.0093589, 1.4870201, -0.1193149, 0.2578262)

imh <- model.normal.nonconjugate.imh(
    phi.0 = c(mean(data), var(data))
  , imh.samples = 100000
  , imh.burnin = 1000
  , alpha = 2
  , beta = 2
  , m = 1
  , v = 10
  , y = data
    )

# Theta/Sigma2 estimates
mean(imh$phi$theta)
mean(imh$phi$sigma2)

# MCMC traces
plot.mcmc.trace(imh$phi$theta, imh$phi$burnin)
plot.mcmc.trace(imh$phi$sigma2, imh$phi$burnin)

# ACF
acf(imh$phi$theta)
acf(imh$phi$sigma2)

# Posterior probability that mu is bigger than 0.5
mean(imh$phi$theta > 0.5)

# MCMC trace on parameter space
plot.mcmc.trace.params(imh$phi$theta, imh$phi$sigma2, imh$phi$burnin)

# Marginal distribution of theta
hist(imh$phi$theta[imh$phi$burnin == "sample"], main = "", xlab = expression(theta), prob = TRUE)
lines(density(imh$phi$theta[imh$phi$burnin == "sample"]), col = "blue")

# Marginal distribution of sigma2
hist(imh$phi$sigma2[imh$phi$burnin == "sample"], main = "", xlab = expression(theta), prob = TRUE)
lines(density(imh$phi$sigma2[imh$phi$burnin == "sample"]), col = "blue")
