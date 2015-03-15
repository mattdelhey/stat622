sample.yij.unconditional <- function(I, J, mu, sigma2.epsilon, sigma2.alpha) {    
    sapply(1:I, function(i) rnorm(J, mu, sqrt(sigma2.alpha + sigma2.epsilon)))
}

sample.yij.conditional <- function(I, J, mu, alpha.i, sigma2.epsilon, sigma2.alpha) {
    alpha.i <- rnorm(I, 0, sqrt(sigma2.alpha))
    sapply(alpha.i, function(a) rnorm(J, mu + a, sqrt(sigma2.epsilon)))
}

model.oneway.anova.theta <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(Y)
    I <- ncol(Y)

    Phi <- matrix(rep(NA, (iters+1)*(I+1)), nrow = iters+1, ncol = I+1)
    Phi[1, ] <- c(init["mu"], rep(init["alpha.i"], I))

    for (i in 2:(iters+1)) {
        mu <- mu.params(ybar.i, sigma2.epsilon, sigma2.alpha)
        mu.new <- rnorm(1, mean = mu[1,1], sd = sqrt(mu[1,2]))
        
        theta <- theta.params(ybar.i, mu.new, sigma2.epsilon, sigma2.alpha)       
        theta.new <- sample.theta(1, theta$theta.j, theta$V)

        Phi[i, ] <- c(mu.new, theta.new)
    }

    return(Phi[-1, ])    
}

model.oneway.anova.eta <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(Y)
    I <- ncol(Y)
    J <- nrow(Y)

    Phi <- matrix(rep(NA, (iters+1)*(I+1)), nrow = iters+1, ncol = I+1)
    Phi[1, ] <- c(init["mu"], rep(init["alpha.i"], I))

    for (i in 2:(iters+1)) {
        eta.new <- sample.eta(1, ybar.i, J, sigma2.epsilon, sigma2.alpha)                
        mu.new <- sample.mu(1, eta = eta.new, sigma2.alpha)
        Phi[i, ] <- c(mu.new, eta.new)
    }

    return(Phi[-1, ])        
}

model.oneway.anova.alpha <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(Y)
    I <- ncol(Y)
    J <- nrow(Y)

    Phi <- matrix(rep(NA, (iters+1)*(I+1)), nrow = iters+1, ncol = I+1)
    Phi[1, ] <- c(init["mu"], rep(init["alpha.i"], I))

    for (i in 2:(iters+1)) {
    }

    return(Phi[-1, ])    
}

sample.mu.alpha <- function(size, ybar.., alphabar, I, J, sigma2.epsilon) {
#    rnorm(size, ybar.. - alphabar, sigma2.epsilon/))
}  

sample.mu <- function(size, eta, sigma2.alpha) {
    J <- length(eta)
    rnorm(size, sum(eta)/J, sqrt(sigma2.alpha/J))
}

sample.eta.j <- function(size, ybar.i, J, sigma2.epsilon, sigma2.alpha) {
    var <- 1 / ( J/sigma2.epsilon + 1/sigma2.alpha )
    mean <- (J*ybar.i / sigma2.epsilon + 1/sigma2.alpha) / ( 1/var )
    rnorm(size, mean, sqrt(var))
}

sample.eta <- function(size, ybar.i, J, sigma2.epsilon, sigma2.alpha) {
    sapply(1:length(ybar.i), function(i)
        sample.eta.j(size, ybar.i[i], J, sigma2.epsilon, sigma2.alpha))
}

mu.params <- function(ybar.i, sigma2.epsilon, sigma2.alpha) {
    V.mu <- 1 / ( sum(1/(sigma2.epsilon/J + sigma2.alpha)) )
    mu.hat <- sum(1/ ( (sigma2.epsilon/J + sigma2.alpha)*ybar.i )) / ( 1/V.mu )
    return(cbind(mu.hat = mu.hat, V.mu = V.mu))
}

theta.params <- function(ybar.i, mu, sigma2.epsilon, sigma2.alpha) {
    J <- length(ybar.i)
    V <- 1 / (J/sigma2.epsilon + 1/sigma2.alpha)
    theta.j <- ( (J/sigma2.epsilon)*ybar.i + (1/sigma2.alpha)*mu ) / ( 1/V )
    return(list(theta.j = theta.j, V = V))
}

sample.theta <- function(size, theta.j, V) {
    sapply(1:length(theta.j), function(j) rnorm(size, theta.j[j], V))
}
