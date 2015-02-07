model.normal.conjuage <- function(x, mu.0, sigma2.0, kappa.0, nu.0) {
    # Summary statistics from data
    y.bar <- mean(x)
    s2 <- var(x)
    n.obs <- length(x)
    ss <- (n.obs - 1) * s2

    # Posterior parameters
    nu.n <- nu.0 + n.obs
    kappa.n <- kappa.0 + n.obs
    mu.n <- (kappa.0*mu.0 + n.obs * y.bar) / kappa.n
    sigma2.n <- (1/nu.n) * (nu.0 * sigma2.0 + (n.obs - 1) *
                              s2 + (kappa.0*n.obs / kappa.n) * (y.bar - mu.0)^2)

    # Calculate SS.n two ways, check for equivalence
    ss.n.2 <- nu.n * sigma2.n
    ss.n <- (nu.0*sigma2.0) + ss + (n.obs*kappa.0)/kappa.n * (y.bar - mu.0)^2
    stopifnot(all.equal(ss.n, ss.n.2))

    posterior.parameters <- list(
        nu.n     = nu.n
      , kappa.n  = kappa.n
      , mu.n     = mu.n
      , sigma2.n = sigma2.n
      , ss.n     = ss.n
      , data     = list(
            y.bar = y.bar
          , s2    = s2
          , n.obs = n.obs
          , ss    = ss))
    return(posterior.parameters)
}

model.normal.conjugate.mc <- function(posterior.parameters, mc.samples) {
    sigma2.posterior.sample <- 1 / rgamma(
        mc.samples
      , shape = posterior.parameters$nu.n / 2
      , rate  = posterior.parameters$sigma2.n * posterior.parameters$nu.n / 2)
    
    theta.posterior.sample <- rnorm(
        mc.samples
      , mean = posterior.parameters$mu.n
      , sd   = sqrt(sigma2.posterior.sample / posterior.parameters$kappa.n))
    
    predictive.posterior.sample <- rnorm(
        mc.samples
      , mean = theta.posterior.sample
      , sd   = sqrt(sigma2.posterior.sample))

    posterior.samples <- list(
        sigma2.posterior.sample     = sigma2.posterior.sample
      , theta.posterior.sample      = theta.posterior.sample
      , predictive.posterior.sample = predictive.posterior.sample)
    return(posterior.samples)
}

mc.quantile.ci <- function(posterior.sample, alpha = 0.05) {
    posterior.mean <- mean(posterior.sample)
    posterior.quantiles <- quantile(posterior.sample, c(alpha/2, 1 - alpha/2))
    
    mc.quantile.ci <- list(
        posterior.mean = posterior.mean
      , posterior.qi.lower = posterior.quantiles[1]
      , posterior.qi.upper = posterior.quantiles[2])
    return(mc.quantile.ci)
}
