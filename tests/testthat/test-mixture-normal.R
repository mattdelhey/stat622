set.seed(1)
glucose <- scan("/Users/mdelhey/Dropbox/Courses/Stat622/hw3/glucose.dat")

alpha <- 1
beta <- 1
mu.0 <- 120
tau2.0 <- 200
sigma2.0 <- 1000
nu.0 <- 10

devtools::load_all("~/stat622")

g <- model.mixture.normal.gibbs(
    phi.0 = c(mean(glucose), mean(glucose), var(glucose), var(glucose), 0.5, length(glucose)/2)
  , gibbs.samples = 10000
  , gibbs.burnin = 1000
  , alpha = alpha
  , beta = beta
  , mu.0 = mu.0
  , tau2.0 = tau2.0
  , nu.0 = nu.0
  , sigma2.0 = sigma2.0
  , y.bar = mean(glucose)
  , n.obs = length(glucose)
  , s2 = var(glucose)
  , y = glucose
    )

str(g)

# MCMC traces
plot.mcmc.trace(g$phi$p, g$phi$burnin)
plot.mcmc.trace(g$phi$sum.x, g$phi$burnin)
