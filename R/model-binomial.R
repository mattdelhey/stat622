model.binomial.conjugate <- function(x = NULL, n.obs = NULL, y.sum = NULL, alpha, beta) {
    # Should I calculate sufficient statistics?
    data <- if (!is.null(x)) model.normal.summary(x) else list(n.obs = n.obs, y.sum = y.sum)

    # Posterior parameters
    alpha.n <- alpha + data$y.sum
    beta.n <- beta + n.obs - data$y.sum

    # Ensoure our parameters are within their support
    stopifnot(alpha.n > 0 & beta.n > 0)

    posterior.parameters <- list(
        alpha.n = alpha.n
      , beta.n  = beta.n
      , size    = data$n.obs)
    return(posterior.parameters)
}

model.binomial.conjugate.mc <- function(posterior.parameters, mc.samples) {
    theta.posterior.sample <- rbeta(
        n      = mc.samples
      , shape1 = posterior.parameters$alpha.n
      , shape2 = posterior.parameters$beta.n)
        
    predictive.posterior.sample <- rbetabinom(
        n      = mc.samples
      , size   = posterior.parameters$size
      , shape1 = posterior.parameters$alpha.n
      , shape2 = posterior.parameters$beta.n)  

    posterior.samples <- list(
        theta.posterior.sample      = theta.posterior.sample
      , predictive.posterior.sample = predictive.posterior.sample)
    
    return(posterior.samples)
}


model.binomial.summary <- function(x) {
    y.sum <- sum(x)
    data <- list(
        y.sum = y.sum
      , n.obs = length(x))
    return(data)
}

rbetabinom <- function(n, size, shape1, shape2) {
    rbinom(n, size = size, prob = rbeta(n, shape1, shape2))
}


