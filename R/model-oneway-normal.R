#' @title shrink.estimate
#' @param lambda  regularization parameter in [0, 1], vector of length K or scalar
#' @param ybar vector of group means
#' @details Lambda = 1 is no pooling (all groups different). Lambda = 0 is complete pooling (all
#' groups same). Values of lambda should probably be > n/sum(n).
shrink.estimate <- function(lambda, ybar, n) {
    grand.mean <- mean(ybar)
    lambda * ybar + (1 - lambda) * grand.mean
}

#' @title emperical.bayes
#' "An Introduction to Emperical Bayes Data Analysis", Casella 1985
emperical.bayes <- function(y, g, sigma2) {
    ybar.j <- tapply(y, g, mean)
    ybar.. <- mean(y)

    #emperical.bayes.est(ybar.j, ybar, sigma2 = 2)
    emperical.bayes.aov(y, g)
}

emperical.bayes.aov <- function(y, g) {
    m <- length(unique(g))
    a <- anova(lm(y ~ g))
    T <- a$"F value"[1]

    ybar.j <- tapply(y, g, mean)
    ybar.. <- mean(y)
    
    (m-3)/(m-1) * 1/T * ybar.. + (1- ((m-3)/(m-1)*(1/T))) *ybar.j
}

#' @title emperical.bayes.est
emperical.bayes.est <- function(ybar.j, ybar, sigma2) {
    ss <- sum( (ybar.j - ybar..)^2 )
    ( ((m-3)*sigma2) / ss ) * ybar.. + (1 - ((m-3)*sigma2) / ss) * ybar.j
}


#' @title emperical.oneway
#' Pg. 118 BDA3
#' The problems of this approach are detailed under "Difficulty with a natural non-Bayesian estimate
#' of the hyperparameters". In summary, we are ignoring our uncertainty about mu and tau and so we
#' will be over-confident in our estimates of theta.
emperical.oneway <- function(y, g, size) {
    ybar.j <- tapply(y, g, mean)
    sigma2.j <- tapply(y, g, var) # assume known sigma2
    ybar.. <- sum(ybar.j / sigma2.j) / sum(1/sigma2.j) # pooled estimate

    # Use ANOVA point estimates for tau2 and mu
    a <- anova(lm(y ~ g))
    ms.b <- a[1, "Mean Sq"]
    ms.w <- a[2, "Mean Sq"]
    tau2 <- (ms.b - ms.w) / length(y) * 100
    mu <- ybar..

    theta.hat <- (ybar.j/sigma2.j + mu/tau2) / (1/sigma2.j + 1/tau2)
    V. <- 1 / (1/sigma2.j + 1/tau2)

    theta <- (sigma2.j / (sigma2.j + tau2))*mu +
      (tau2 / (tau2 + sigma2.j)) * ybar.j

    return(theta)
}

#' @title oneway.uniform
#' Pg. 115 BDA3
#' Assumes that sigma2 is known. In this case, use sample estimate.
#' This isn't a horrible assumption, see BDA3 end of Ch 8.
oneway.uniform <- function(y, g, size) {
    ybar.j <- tapply(y, g, mean)
    sigma2.j <- tapply(y, g, var) # assume known sigma2

    tau2 <- sample.oneway.tau2(size, sigma2.max = 120, sigma2.j, ybar.j)
    mu <- sample.oneway.mu(size, tau2, ybar.j, sigma2.j)
    theta.j <- sample.oneway.theta(size, mu, tau2, ybar.j, sigma2.j)
    
    colnames(theta.j) <- names(ybar.j)
    return(list(tau2 = tau2, mu = mu, theta.j = theta.j))
}

#' @title sample.oneway.theta
sample.oneway.theta <- function(size, mu, tau2, ybar.j, sigma2.j) {
    theta.est <- t(sapply(1:size, function(i)
        oneway.theta(mu[i], tau2[i], ybar.j, sigma2.j)))
            
    t(sapply(1:size, function(i)
        sapply(1:length(ybar.j), function(j)
            rnorm(1, unlist(theta.est[i, 1])[j], sqrt(unlist(theta.est[i, 2])[j])))))
}

#' @title sample.oneway.mu
sample.oneway.mu <- function(size, tau2, ybar.j, sigma2.j) {
    mu.est <- t(sapply(tau2, oneway.mu, ybar.j, sigma2.j))
    sapply(1:size, function(i)
        rnorm(1, mu.est[i, 1], sqrt(mu.est[i, 2])))
}

#' @title sample.oneway.tau2
sample.oneway.tau2 <- function(size, sigma2.max, sigma2.j, ybar.j) {
    tau2.grid <- seq(from = 1, to = sigma2.max, length.out = 1000)
    mu.est <- t(sapply(tau2.grid, oneway.mu, ybar.j, sigma2.j))

    p.tau2 <- sapply(1:length(tau2.grid), function(i)
        p.tau.given.y(tau2 = tau2.grid[i], mu.est = mu.est[i, 1], V.mu = mu.est[i, 2],
                      sigma2.j = sigma2.j, ybar.j = ybar.j))
    #p.tau2 <- sapply(tau2.grid, p.tau.given.y, mu.est[, 1], mu.est[, 2], sigma2.j, ybar.j)
    #p.tau2 <- sapply(tau2.grid, oneway.tau2, sigma2.j, ybar.j)
    sample(tau2.grid, size, prob = exp(p.tau2 - max(p.tau2)), replace = TRUE)
}

#' @title p.tau.given.y
#' @note Should be changed to log-scale for numerical stability.
p.tau.given.y <- function(tau2, mu.est, V.mu, ybar.j, sigma2.j, p.tau = 1 ) {
    #log(p.tau) + 0.5*log(V.mu) +
    #  sum( -0.5*log(sigma2.j + tau2) + (-(ybar.j - mu.est)^2 / (2*(sigma2.j +tau2))) )
    sum(dnorm(ybar.j, mu.est, sqrt(sigma2.j + tau2), log = TRUE)) -
      dnorm(mu.est, mu.est, sqrt(V.mu), log = TRUE)
}
  
#' @title oneway.tau2
#' Pg. 117 BDA3
#' @note Assumes that tau2 \propto 1
oneway.tau2 <- function(tau2, sigma2.j, ybar.j) {
    #est <- oneway.mu(tau2, ybar.j, sigma2.j, list = TRUE)
    est <- t(sapply(tau2, oneway.mu, ybar.j, sigma2.j))
    #sqrt(est$V.mu) * prod( (sigma2.j + tau2)^(-1/2) *
    #  exp( -(ybar.j - est$mu.hat)^2 / (2*sigma2.j + 2*tau2) ) )
    sqrt(est[, 2]) * prod( (sigma2.j + tau2)^(-1/2) *
      exp( -(ybar.j - est[, 1])^2 / (2*sigma2.j + 2*tau2) ) )
}

#' @title oneway.theta
#' Pg. 116 BDA3
oneway.theta <- function(mu, tau2, ybar.j, sigma2.j) {
    theta.hat.j <- (ybar.j/sigma2.j + mu/tau2) / (1/sigma2.j + 1/tau2)
    V.j <- 1 / (1/sigma2.j + 1/tau2)
    #theta.hat.j <- outer(mu/tau2, ybar.j/sigma2.j, "+") / outer(1/tau2, 1/sigma2.j, "+")
    #V.j <- t(1 / outer(1/sigma2.j, 1/tau2, "+"))
    return(list(theta.hat.j = theta.hat.j, V.j = V.j))
}

#' @title oneway.mu
#' Pg. 117 BDA3
oneway.mu <- function(tau2, ybar.j, sigma2.j) {
    mu.hat <- sum(ybar.j / (sigma2.j + tau2)) / sum(1 / (sigma2.j + tau2))
    V.mu <- 1 / sum(1 / (sigma2.j + tau2))
    return(cbind(mu.hat, V.mu))        
    ## mu.hat <- colSums(matrix(1/ybar.j, ncol = 25) %*% outer(sigma2.j, tau2, "+")) /
    ##   colSums(1 / outer(sigma2.j, tau2, "+"))
    ## V.mu <- 1 / colSums(1 / outer(sigma2.j, tau2, "+"))
    ##return(list(mu.hat = mu.hat, V.mu = V.mu))
}

#' @title oneway.normal.gibbs.sampler
#' Gibbs sampler for the "one-way normal random effects model with known
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
oneway.gibbs.sampler <- function(y = NULL, g, n.samples, n.burnin, priors,
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

mowd <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
