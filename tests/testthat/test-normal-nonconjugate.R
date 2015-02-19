set.seed(1)

alpha <- 2
beta  <- 2
m     <- 1
v     <- 10

imh.samples <- 10000
imh.burnin  <- 1000

data <- c(2.3656491, 2.4952035, 1.0837817, 0.7586751, 0.8780483, 1.2765341,
          1.4598699, 0.1801679, -1.0093589, 1.4870201, -0.1193149, 0.2578262)

imh <- model.normal.nonconjugate.imh(
    phi.0 = c(mean(data), var(data))
  , imh.samples = imh.samples
  , imh.burnin  = imh.burnin
  , alpha       = alpha
  , beta        = beta
  , m           = m
  , v           = v
  , y           = data
    )

# Theta/Sigma2 estimates
mean(imh$phi$theta)
mean(imh$phi$sigma2)

# test plotting of MCMC traces
p <- plot.mcmc.trace(imh$phi$theta, imh$phi$burnin)
p <- plot.mcmc.trace(imh$phi$sigma2, imh$phi$burnin)

# test ACF
a <- acf(imh$phi$theta)
a <- acf(imh$phi$sigma2)

# test MCMC trace on parameter space
p <- plot.mcmc.trace.params(imh$phi$theta, imh$phi$sigma2, imh$phi$burnin)

# test marginal posterior distribution plots
plot.mcmc.marginal(imh$phi$theta, imh$phi$burnin, main = "marginal posterior of theta")
plot.mcmc.marginal(imh$phi$sigma2, imh$phi$burnin, main = "marginal posterior of sigma2", xlab = "sigma2")
