model.mixture.beta <- function(x, alpha1) {
    data <- model.normal.summary(x)
    return(posterior.parameters)
}

model.mixture.beta.importance <- function(n.samples, h, alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2, ...) {
    if (n.samples %% 1 != 0) stop("number of samples must be an integer")
    if (!is.function(h)) stop("h must be a function (use identity for expectation)")

    #x <- rbeta(n.samples, 1, 1)
    #weights <- dbetamixture(x, alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2) / dbeta(x, 1, 1)

    x <- runif(n.samples, 0, 1)
    weights <- dbetamixture(x, alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2) / dunif(x, 0, 1)
    
    weights.std <- weights / sum(weights)
    expectation.h <- sum(h(x, ...) * weights.std)
    return(expectation.h)
}

rbetamixture <- function(n, alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2) {
    check.betamixture(alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2)
    delta.1*rbeta(n, alpha.1, beta.1) + delta.2*rbeta(n, alpha.2, beta.2)   
}

dbetamixture <- function(x, alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2) {
    check.betamixture(alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2)
    delta.1*dbeta(x, alpha.1, beta.1) + delta.2*dbeta(x, alpha.2, beta.2)   
}

check.betamixture <- function(alpha.1, beta.1, delta.1, alpha.2, beta.2, delta.2) {
    stopifnot(delta.1 + delta.2 == 1 & alpha.1 > 0 & alpha.2 > 0 & beta.1 > 0 & beta.2 > 0)
}


    
