#' @title model.mixture
#' Posterior parameters for the semiconjugate mixture model
model.mixture <- function(x, ) {
    
}

model.mixture.gibbs <- function(phi.0, gibbs.iters,
                                a, b, mu.0, tau2.0, nu.0, sigma2.0,
                                y.bar, n.obs, s2) {
    if (length(phi.0) !=  & !is.vector(phi.0)) stop("phi must be a vector of length 2")
    if (gibbs.iters %% 1 != 0) stop("number of iterations must be an integer")

    # Initialize phi [matrix of dependent sequence of posterior samples]
    phi.empty <- data.frame(
        theta.1 = rep(NA, gibbs.iters)
      , theta.2 = rep(NA, gibbs.iters)
      , sigma2.1 = rep(NA, gibbs.iters)
      , sigma2.2 = rep(NA, gibbs.iters)
      , p = rep(NA, gibbs.iters)
      , sum.x = rep(NA, gibbs.iters))
    phi <- rbind(phi.0, phi.empty)

    for (i in 2:(gibbs.iters+1)) {
        # Full conditional of p
        p.conditional <- model.binomial.conjugate(
            y.sum = phi$sumx[i-1]
          , alpha = a
          , beta = b
          , n.obs = n.obs)
        phi$p[i] <- rbeta(1, p.conditional$alpha.n, p.conditional$beta.n)

        # Full conditional of X
        phi$sum.x[i] <- rbinom(1, n.obs, phi$p[i] * rnorm(1, mean = phi$theta.1[i-1], sd = sqrt(phi$sigma.1[i-1])))

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
    return(phi)
}


normal.mixture.conjugate <- function()
