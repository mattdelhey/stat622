devtools::load_all("~/stat622")
set.seed(1)
glucose <- scan("/Users/mdelhey/Dropbox/Courses/Stat622/hw3/glucose.dat")

alpha <- 1
beta <- 1
mu.0 <- 120
tau2.0 <- 200
sigma2.0 <- 1000
nu.0 <- 10

gibbs.samples <- 1000
gibbs.burnin <- 100

system.time(
g <- model.mixture.normal.gibbs(
    phi.0 = c(
        theta.1  = mean(glucose)
      , theta.2  = mean(glucose)
      , sigma2.1 = var(glucose)
      , sigma2.2 = var(glucose)
      , p        = 0.5
      , sum.x    = length(glucose)/2)
  , gibbs.samples = gibbs.samples
  , gibbs.burnin  = gibbs.burnin
  , alpha    = alpha
  , beta     = beta
  , mu.0     = mu.0
  , tau2.0   = tau2.0
  , nu.0     = nu.0
  , sigma2.0 = sigma2.0
  , y        = glucose
  , verbose  = TRUE
    )
    )

str(g)

mean(g$phi$sum.x)
mean(rowSums(g$x.mat))

# MCMC traces
plot.mcmc.trace(g$phi$p, g$phi$burnin)
plot.mcmc.trace(g$phi$sum.x, g$phi$burnin)
plot.mcmc.trace(g$phi$theta.1, g$phi$burnin)
plot.mcmc.trace(g$phi$sigma2.2, g$phi$burnin)

# Marginal of sum(x)
plot.mcmc.marginal(g$phi$sum.x, g$phi$burnin)

# Expectations
mean(g$phi$p)
mean(g$phi$theta.1)
mean(g$phi$theta.2)
mean(g$phi$sigma2.1)
mean(g$phi$sigma2.2)

