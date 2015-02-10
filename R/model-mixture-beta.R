model.mixture.beta.importance <- function(n.samples, h, alpha.1, beta.1, alpha.2, beta.2, delta, ...) {
    if (!is.function(h)) stop("h must be a function (use identity for expectation)")

    #x <- rbeta(n.samples, 1, 1)
    #weights <- dbetamixture(x, alpha.1, beta.1, alpha.2, beta.2, delta) / dbeta(x, 1, 1)

    x <- runif(n.samples, 0, 1)
    weights <- dbetamixture(x, alpha.1, beta.1, alpha.2, beta.2, delta) / dunif(x, 0, 1)

    # Standardized weights
    weights.std <- weights / sum(weights)
    expectation.h <- sum(h(x, ...) * weights.std)    

    # Raw weights
    expectation.h.raw <- sum(h(x, ...) * weights) / n.samples

    if (!isTRUE(all.equal(expectation.h, expectation.h.raw, tolerance = 0.01)))
        warning("Expectation differs for standardized and raw weights ")
    
    return(expectation.h)
}

rbetamixture <- function(n, alpha.1, beta.1, alpha.2, beta.2, delta) {
    check.betamixture(n, alpha.1, beta.1, alpha.2, beta.2, delta)
    state.vec <- sample(1:2, prob = c(delta, 1 - delta), size = n.samples, replace = TRUE)    
    rbeta(n = n.samples, shape1 = c(alpha.1, alpha.2)[state.vec], shape2 = c(beta.1, beta.2)[state.vec])
}

rbetamixtrue

dbetamixture <- function(x, alpha.1, beta.1, alpha.2, beta.2, delta) {
    check.betamixture(length(x), alpha.1, beta.1, alpha.2, beta.2, delta)
    delta*dbeta(x, alpha.1, beta.1) + (1-delta)*dbeta(x, alpha.2, beta.2)   
}

check.betamixture <- function(n, alpha.1, beta.1, alpha.2, beta.2, delta) {
    stopifnot(delta >= 0 & delta <= 1 & alpha.1 > 0 & alpha.2 > 0 & beta.1 > 0 & beta.2 > 0 & n %% 1 == 0)
}


    
