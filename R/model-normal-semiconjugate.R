#' @title model.normal.semiconjugate.gibbs
#' Gibbs samples for the semiconjugate normal model
model.normal.semiconjugate.gibbs <- function(phi.0, gibbs.samples, gibbs.burnin = 0,
                                             mu.0, tau2.0, nu.0, sigma2.0,
                                             y.bar, n.obs, s2) {
    if (gibbs.samples %% 1 != 0 | gibbs.burnin %% 1 != 0)
        stop("number of samples must be an integer")

    gibbs.iters <- gibbs.samples + gibbs.burnin + 1
    phi <- construct.phi(phi0, gibbs.iters, gibbs.burnin, vars = c("theta", "sigma2"))

    # Account for phi.0 by starting at 2
    for (i in 2:gibbs.iters) {
        # Use new sigma2 to calculate theta posterior parameters
        theta.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2[i-1]
          , mu.0 = mu.0
          , tau2.0 = tau2.0
          , n.obs = n.obs
          , y.bar = y.bar)

        # Generate new theta
        phi$theta[i] <- rnorm(1, theta.conditional$mu.n,
                              sqrt(theta.conditional$tau2.n))

        # Use new theta to calculate sigma2 posterior parameters
        sigma2.conditional <- model.normal.semiconjugate.sigma2(
            theta = phi$theta[i]
          , nu.0 = nu.0
          , sigma2.0 = sigma2.0
          , n.obs = n.obs
          , y.bar = y.bar
          , s2 = s2)
        
        # Generate new sigma2
        phi$sigma2[i] <- 1 / rgamma(1, sigma2.conditional$nu.n / 2,
                                sigma2.conditional$nu.n * sigma2.conditional$sigma2.n / 2)
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
      , tau2.n = tau2.n)
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
      , sigma2.n = sigma2.n)
    return(posterior.parameters.sigma2)
}
