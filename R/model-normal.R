#' @title model.normal.conjugate
#' Posterior parameters for the conjugate normal model
model.normal.conjuage <- function(x, mu.0, sigma2.0, kappa.0, nu.0) {
    data <- model.normal.summary(x)

    # Posterior parameters
    nu.n <- nu.0 + data$n.obs
    kappa.n <- kappa.0 + data$n.obs
    mu.n <- (kappa.0*mu.0 + data$n.obs * data$y.bar) / kappa.n
    sigma2.n <- (1/nu.n) * (nu.0 * sigma2.0 + (data$n.obs - 1) *
                              data$s2 + (kappa.0*data$n.obs / kappa.n) * (data$y.bar - mu.0)^2)

    # Calculate SS.n two ways, check for equivalence
    ss.n.2 <- nu.n * sigma2.n
    ss.n <- (nu.0*sigma2.0) + ss + (n.obs*kappa.0)/kappa.n * (y.bar - mu.0)^2
    stopifnot(all.equal(ss.n, ss.n.2))

    posterior.parameters <- list(
        nu.n     = nu.n
      , kappa.n  = kappa.n
      , mu.n     = mu.n
      , sigma2.n = sigma2.n
      , ss.n     = ss.n)
    return(posterior.parameters)
}

#' @title model.normal.conjugate.mc
#' Monte Carlo samples for the conjugate normal model
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

#' @title mc.quantile.ci
#' Monte Carlo quantile-based confidence interval
mc.quantile.ci <- function(posterior.sample, alpha = 0.05) {
    posterior.mean <- mean(posterior.sample)
    posterior.quantiles <- quantile(posterior.sample, c(alpha/2, 1 - alpha/2))
    
    mc.quantile.ci <- list(
        posterior.mean = posterior.mean
      , posterior.qi.lower = posterior.quantiles[1]
      , posterior.qi.upper = posterior.quantiles[2])
    return(mc.quantile.ci)
}

#' @title model.normal.summary
#' Summary statistics for normal data
model.normal.summary <- function(x) {
    y.bar <- mean(x)
    s2 <- var(x)
    n.obs <- length(x)
    ss <- (n.obs - 1) * s2

    data <- list(
        y.bar = y.bar
      , s2    = s2
      , n.obs = n.obs
      , ss    = ss)
    return(data)
}

#' @title model.normal.semiconjugate.gibbs
#' Gibbs samples for the semiconjugate normal model
model.normal.semiconjugate.gibbs <- function(phi.0, gibbs.iters,
                                             mu.0, tau2.0, nu.0, sigma2.0,
                                             y.bar, n.obs, s2) {
    if (length(phi.0) != 2 & !is.vector(phi.0)) stop("phi must be a vector of length 2")
    if (gibbs.iters %% 1 != 0) stop("number of iterations must be an integer")

    # Initialize phi [matrix of dependent sequence of posterior samples]
    phi.empty <- data.frame(theta = rep(NA, gibbs.iters), sigma2 = rep(NA, gibbs.iters))
    phi <- rbind(phi.0, phi.empty)

    for (i in 2:(gibbs.iters+1)) {
        # Use new sigma2 to calculate theta posterior parameters
        theta.posterior <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2[i-1]
          , mu.0 = mu.0
          , tau2.0 = tau2.0
          , n.obs = n.obs
          , y.bar = y.bar
            )

        # Generate new theta
        phi$theta[i] <- rnorm(1, theta.posterior$mu.n,
                              sqrt(theta.posterior$tau2.n))

        # Use new theta to calculate sigma2 posterior parameters
        sigma2.posterior <- model.normal.semiconjugate.sigma2(
            theta = phi$theta[i]
          , nu.0 = nu.0
          , sigma2.0 = sigma2.0
          , n.obs = n.obs
          , y.bar = y.bar
          , s2 = s2
            )

        # Generate new sigma2
        phi$sigma2[i] <- 1 / rgamma(1, sigma2.posterior$nu.n / 2,
                                sigma2.posterior$nu.n * sigma2.posterior$sigma2.n / 2)
    }    
    return(phi)
}

#' @title model.normal.semiconjugate.theta
#' Full conditional distribution for theta (sigma2 known)
model.normal.semiconjugate.theta <- function(sigma2, mu.0, tau2.0, n.obs, y.bar) {
    mu.n <- ( mu.0/tau2.0 + n.obs*y.bar/sigma2 ) / ( 1/tau2.0 + n.obs/sigma2 )
    tau2.n <- 1 / ( 1/tau2.0 + n.obs/sigma2 )
    posterior.parameters.theta <- list(
        mu.n = mu.n
      , tau2.n = tau2.n
        )
    return(posterior.parameters.theta)
}

#' @title model.normal.semiconjugate.sigma2
#' Full conditional distribution for sigma2 (theta known)
model.normal.semiconjugate.sigma2 <- function(theta, nu.0, sigma2.0, n.obs, y.bar, s2) {
    nu.n <- nu.0 + n.obs
    # Use trick to calculate sample variance from summary statistics
    sigma2.n <- (1/nu.n) * (nu.0*sigma2.0 + (n.obs-1)*s2 + n.obs*(y.bar - theta)^2)
    posterior.parameters.sigma2 <- list(
        nu.n = nu.n
      , sigma2.n = sigma2.n
        )
    return(posterior.parameters.sigma2)
}
