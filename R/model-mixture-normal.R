model.mixture.normal.gibbs <- function(phi.0, gibbs.samples, gibbs.burnin,
                                       alpha, beta, mu.0, tau2.0, nu.0, sigma2.0,
                                       y, verbose = FALSE) {
    if (gibbs.samples %% 1 != 0) stop("number of iterations must be an integer")

    n.obs <- length(y)

    gibbs.iters <- gibbs.samples + gibbs.burnin + 1
    phi <- construct.phi(phi.0, gibbs.iters, gibbs.burnin, vars = names(phi.0))

    for (i in 2:gibbs.iters) {
        if (verbose & i %% 500 == 0)
            message(sprintf("iteration: %i", i))
        
        # Full conditional (posterior) of p
        p.conditional <- model.binomial.conjugate(
            y.sum = phi$sum.x2[i-1]
          , alpha = alpha
          , beta  = beta
          , n.obs = n.obs)
        phi$p[i] <- rbeta(1, p.conditional$alpha.n, p.conditional$beta.n)
        
        # Full conditional of X
        z.1 <- prod(dnorm(y, mean = phi$theta.1[i-1], sd = sqrt(phi$sigma2.1[i-1]), log = FALSE))
        z.2 <- prod(dnorm(y, mean = phi$theta.2[i-1], sd = sqrt(phi$sigma2.2[i-1]), log = FALSE))
        phi$sum.x[i] <- rbinom(1, n.obs, phi$p[i]*z.1 / (phi$p[i]*z.1 + (1-phi$p[i])*z.2))

        #x <- mclapply(1:n.obs, function(j) {
        x <- parLapply(cl, 1:n.obs, function(j) {
            dn1 <- dnorm(y[j], mean = phi$theta.1[i-1], sd = sqrt(phi$sigma2.1[i-1]))
            dn2 <- dnorm(y[j], mean = phi$theta.2[i-1], sd = sqrt(phi$sigma2.2[i-1]))
            rbinom(1, 1, phi$p[i]*dn1 / (phi$p[i]*dn1 + (1-phi$p[i])*dn2))
        })
        x.vec <- unlist(x)
        phi$sum.x2[i] <- sum(x.vec)

        # Define new data statistics for Y_1 and Y_2
        n.obs.1 <- sum(x.vec)
        n.obs.2 <- n.obs - n.obs.1
        y.bar.1 <- sum(y[x.vec == 1]) / n.obs.1
        y.bar.2 <- sum(y[x.vec == 0]) / n.obs.2
        s2.1 <- var(y[x.vec == 1])
        s2.2 <- var(y[x.vec == 0])
        stopifnot(n.obs.1 + n.obs.2 == n.obs)

        # Full conditional of theta.1
        theta.1.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2.1[i-1]
          , mu.0   = mu.0
          , tau2.0 = tau2.0
          , n.obs  = n.obs.1
          , y.bar  = y.bar.1
            )
        phi$theta.1[i] <- rnorm(1, theta.1.conditional$mu.n,
                                sqrt(theta.1.conditional$tau2.n))

        # Full conditional of theta.2
        theta.2.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2.2[i-1]
          , mu.0   = mu.0
          , tau2.0 = tau2.0
          , n.obs  = n.obs.2
          , y.bar  = y.bar.2
            )
        phi$theta.2[i] <- rnorm(1, theta.2.conditional$mu.n,
                                sqrt(theta.2.conditional$tau2.n))

        # Use new theta to calculate sigma2 posterior parameters
        sigma2.1.conditional <- model.normal.semiconjugate.sigma2(
            theta    = phi$theta.1[i]
          , nu.0     = nu.0
          , sigma2.0 = sigma2.0
          , n.obs    = n.obs.1
          , y.bar    = y.bar.1
          , s2       = s2.1
            )
        phi$sigma2.1[i] <- 1 / rgamma(1, sigma2.1.conditional$nu.n / 2,
                                      sigma2.1.conditional$nu.n * sigma2.1.conditional$sigma2.n / 2)
        
        sigma2.2.conditional <- model.normal.semiconjugate.sigma2(
            theta    = phi$theta.2[i]
          , nu.0     = nu.0
          , sigma2.0 = sigma2.0
          , n.obs    = n.obs.2
          , y.bar    = y.bar.2
          , s2       = s2.2
            )
        phi$sigma2.2[i] <- 1 / rgamma(1, sigma2.2.conditional$nu.n / 2,
                                      sigma2.2.conditional$nu.n * sigma2.2.conditional$sigma2.n / 2)
    }
    
    phi <- list(
        phi = phi
      , eff.theta.1  = as.numeric(coda::effectiveSize(phi$theta.1))
      , eff.theta.2  = as.numeric(coda::effectiveSize(phi$theta.2))
      , eff.sigma2.1 = as.numeric(coda::effectiveSize(phi$sigma2.1))
      , eff.sigma2.2 = as.numeric(coda::effectiveSize(phi$sigma2.1))
      , eff.p        = as.numeric(coda::effectiveSize(phi$p))
      , eff.sum.x    = as.numeric(coda::effectiveSize(phi$sum.x))
      , eff.sum.x2   = as.numeric(coda::effectiveSize(phi$sum.x2)))
      #, eff2.sum.x2  = effective.sample.size(phi$sum.x2, n.lags = 1000))
    return(phi)
}


rnormalmixture <- function(n, theta.1, theta.2, sigma2.1, sigma2.2, delta) {
    check.normalmixture(n, theta.1, theta.2, sigma2.1, sigma2.2, delta)
    # Equivilent to binomial(1, delta)
    state.vec <- sample(1:2, prob = c(delta, 1 - delta), size = n.samples, replace = TRUE)
    # State vectors
    theta.vec <- c(theta.1, theta.2)
    sigma2.vec <- c(sigma2.1, sigma2.2)
    # Sample from normal with delta-proportions
    rnorm(n = n, mu = theta.vec[state.vec], sd = sqrt(sigma2.vec[state.vec]))
}

dbetamixture <- function(x, theta.1, theta.2, sigma2.1, sigma2.2, delta) {
    check.normalmixture(length(x), theta.1, theta.2, sigma2.1, sigma2.2, delta)
    delta*dnorm(x, theta.1, sqrt(sigma2.1)) + (1-delta)*dnorm(x, theta.2, sqrt(sigma2.2))
}

check.normalmixture <- function(n, theta.1, theta.2, sigma2.1, sigma2.2, delta) {
    stopifnot(n %% 1 == 0, sigma2.1 > 0, sigma2.2 > 0, is.numeric(theta.1), is.numeric(theta.2),
              delta >= 0, delta <= 1)
}
