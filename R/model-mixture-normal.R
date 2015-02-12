model.mixture.normal.gibbs <- function(phi.0, gibbs.samples, gibbs.burnin,
                                       alpha, beta, mu.0, tau2.0, nu.0, sigma2.0,
                                       y.bar, n.obs, s2, y) {
    if (gibbs.samples %% 1 != 0) stop("number of iterations must be an integer")

    gibbs.iters <- gibbs.samples + gibbs.burnin + 1
    phi <- construct.phi(phi.0, gibbs.iters, gibbs.burnin,
                         vars = c("theta.1", "theta.2", "sigma2.1", "sigma2.2", "p", "sum.x"))

    for (i in 2:gibbs.iters) {
        # Full conditional (posterior) of p
        p.conditional <- model.binomial.conjugate(
            y.sum = phi$sum.x[i-1]
          , alpha = alpha
          , beta = beta
          , n.obs = n.obs)
        phi$p[i] <- rbeta(1, p.conditional$alpha.n, p.conditional$beta.n)

        # Full conditional of X
        phi$sum.x[i] <- rbinom(1, n.obs, phi$p[i] *
                                 prod(dnorm(y, mean = phi$theta.1[i-1], sd = sqrt(phi$sigma2.1[i-1]))))

        # Full conditional of theta.1
        theta.1.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2.1[i-1]
          , mu.0 = mu.0
          , tau2.0 = tau2.0
          , n.obs = n.obs
          , y.bar = y.bar
            )
        phi$theta.1[i] <- rnorm(1, theta.1.conditional$mu.n,
                                sqrt(theta.1.conditional$tau2.n))

        # Full conditional of theta.2
        theta.2.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2.2[i-1]
          , mu.0 = mu.0
          , tau2.0 = tau2.0
          , n.obs = n.obs
          , y.bar = y.bar
            )
        phi$theta.2[i] <- rnorm(1, theta.2.conditional$mu.n,
                                sqrt(theta.2.conditional$tau2.n))
        

        # Use new theta to calculate sigma2 posterior parameters
        sigma2.1.conditional <- model.normal.semiconjugate.sigma2(
            theta = phi$theta.1[i]
          , nu.0 = nu.0
          , sigma2.0 = sigma2.0
          , n.obs = n.obs
          , y.bar = y.bar
          , s2 = s2
            )
        phi$sigma2.1[i] <- 1 / rgamma(1, sigma2.1.conditional$nu.n / 2,
                                      sigma2.1.conditional$nu.n * sigma2.1.conditional$sigma2.n / 2)
        
        sigma2.2.conditional <- model.normal.semiconjugate.sigma2(
            theta = phi$theta.2[i]
          , nu.0 = nu.0
          , sigma2.0 = sigma2.0
          , n.obs = n.obs
          , y.bar = y.bar
          , s2 = s2
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
      , eff.sum.x    = as.numeric(coda::effectiveSize(phi$sum.x)))
    return(phi)
}

