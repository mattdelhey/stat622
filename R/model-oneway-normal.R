#' @title oneway.normal.gibbs.sampler
#' @description Gibbs sampler for the "one-way normal random effects model with known
#' variance". This is a simple case of hierarchical linear models.
#' @details
#' The following is the code for the Gibbs sampler. The result is iid samples from the posterior
#' distribution (our beliefs about network lift and its uncertainty after seeing the data).
#' 
#' Its actually a pretty simple idea: write down the full conditional distributions of the parameters
#' and iteraively generate draws from each using their most recent estimates. By certain Markov
#' properties, the draws converage to draws of the posterior (stationary) as $n \rightarrow \infty$.
#'    
#' The model is similar to the Frequentist one-way random effects model (`nlme` and `lme` in R) and
#' assumes that the within-network variance is equal across all networks (which is likely violated in
#' the case of networks like ESPN). 
one.way.gibbs.sampler <- function(y = NULL, g, n.samples, n.burnin, priors,
                                  y.bar = NULL, y.s2 = NULL, y.n = NULL) {

    stopifnot(is.vector(priors), is.vector(g), n.samples %% 1 == 0, n.burnin %% 1 ==0)
    if (is.null(y) && is.null(y.bar))
        stop("you must specify either data or summary statistics")

    # Keep an extra iter for initial values
    n.iters <- n.samples + n.burnin + 1

    groups <- unique(g)     # group names
    m <- length(groups)     # number of groups

    # Data or sufficient statistics?
    use.data <- (is.null(y.bar) || is.null(y.s2) || is.null(y.n))
    n <- if (use.data) length(y) else sum(y.n) # number of obs
       
    # Initalize prior parameters using sample estimates
    if (use.data) { # Estimate from data
        message("Usinig data to calculate summary statistics.")
        y.n <- y.s2 <- y.bar <- rep(NA, m)
        for (j in 1:m) {
            group.data <- y[g == groups[j]]
            y.bar[j]   <- mean(group.data)
            y.s2[j]    <- var(group.data)
            y.n[j]     <- length(group.data) 
        }
    }
    
    # First (init) row of parameter matrix
    phi.0 <- c(y.bar, mean(y.s2), mean(y.bar), var(y.bar))
    theta.names <- if (!is.null(names(g))) names(g)
                   else if (is.character(groups)) groups
                   else sprintf("theta%i", 1:m) 
    phi.names <- c(theta.names, "sigma2", "mu", "tau2")

    # Full parameter matrix [(m + 3) x n.iters]
    # Extra three comes from sigma2, mu, and tau2 which don't vary by group
    phi <- construct.phi(phi.0, n.iters, n.burnin, phi.names)

    for (i in 2:n.iters) {       
        # Sample theta's from full conditionals
        for (j in 1:m) {
            theta.var <- 1 / ( y.n[j]/phi$sigma2[i-1] + 1/phi$tau2[i-1] )
            theta.mean <- theta.var * ( y.bar[j] * y.n[j] / phi$sigma2[i-1] + phi$mu[i-1]/phi$tau2[i-1] )
            phi[i, j] <- rnorm(1, theta.mean, sqrt(theta.var))
        }

        # Sample sigma2 from full conditional
        nu.n <- priors["nu.0"] + n
        ss <- priors["nu.0"] * priors["sigma2.0"]
        for (j in 1:m) {
            ss <- ss + sum( (y[g == groups[j]] - phi[i, j])^2 ) 
        }
        phi$sigma2[i] <- 1 / rgamma(1, nu.n/2, ss/2)

        # Sample mu from full conditional
        mu.var <- 1 / ( m / phi$tau2[i-1] + 1/priors["gamma2.0"] )
        mu.mean <- mu.var * (m * rowMeans(phi[i, 1:m]) / phi$tau2[i-1] + priors["mu.0"] / priors["gamma2.0"] )
        phi$mu[i] <- rnorm(1, mu.mean, sqrt(mu.var))

        # Sample tau2 from full conditional
        eta.m <- priors["eta.0"] + m
        ss <- priors["eta.0"] * priors["tau2.0"] + sum( (phi[i, 1:m] - phi$mu[i])^2 )
        phi$tau2[i] <- 1 / rgamma(1, eta.m/2, ss/2)
    }    
    return(phi)
}

var.ss <- function(sum.x, sum.x.sq, n) {
    if (any(n < 2)) stop("n < 2")   
    var <- (sum.x.sq - (sum.x)^2/n) / (n-1)
    if (any(var < 0)) stop("negative variance")
    var
}

mean.ss <- function(sum.x, n) sum.x / n


Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
