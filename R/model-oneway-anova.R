sample.yij.unconditional <- function(I, J, mu, sigma2.epsilon, sigma2.alpha) {    
    sapply(1:I, function(i) rnorm(J, mu, sqrt(sigma2.alpha + sigma2.epsilon)))
}

sample.yij.conditional <- function(I, J, mu, alpha.i, sigma2.epsilon, sigma2.alpha) {
    alpha.i <- rnorm(I, 0, sqrt(sigma2.alpha))
    sapply(alpha.i, function(a) rnorm(J, mu + a, sqrt(sigma2.epsilon)))
}

model.oneway.anova.alpha <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(Y)
    ybar.. <- mean(ybar.i)
    I <- ncol(Y)
    J <- nrow(Y)
    Phi <- matrix(rep(NA, (iters+1)*(I+1)), nrow = iters+1, ncol = I+1)
    Phi[1, ] <- c(init["mu"], rep(init["alpha.i"], I))    
    for (i in 2:(iters+1)) {
        alpha.new <- sample.alpha(1, ybar.i, mu = Phi[i-1,1], I, J, sigma2.epsilon, sigma2.alpha)
        mu.new <- sample.mu.alpha(1, ybar.., alphabar = mean(alpha.new), I, J, sigma2.epsilon)
        Phi[i, ] <- c(mu.new, alpha.new)
    }
    return(Phi[-1, ])    
}

model.oneway.anova.eta <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(as.matrix(Y))
    ybar.. <- mean(ybar.i)
    I <- ncol(Y)
    J <- nrow(Y)
    Phi <- matrix(rep(NA, (iters+1)*(I)), nrow = iters+1, ncol = I)
    Phi[1, ] <- c(rep(init["alpha.i"], I))
    for (i in 2:(iters+1)) {
        #eta.new <- sample.eta(1, Y, ybar.i, ybar.., I, J, sigma2.epsilon, sigma2.alpha)
        eta.new <- sapply(1:I, function(i) sample.eta.i(1, ybar.i[i], sigma2.alpha))
        #mu.new <- sample.mu(1, eta.new, I, J, sigma2.alpha)
        Phi[i, ] <- c(eta.new)
    }
    return(Phi[-1, ])        
}

sample.alpha <- function(size, ybar.i, mu, I, J, sigma2.epsilon, sigma2.alpha) {
    sapply(1:I, function(i) sample.alpha.i(1, ybar.i[i], mu, I, J, sigma2.epsilon, sigma2.alpha))
}

sample.alpha.i <- function(size, ybar.i, mu, I, J, sigma2.epsilon, sigma2.alpha) {
    var <- 1 / (J/sigma2.epsilon + 1/sigma2.alpha)
    #mean <- ( J*sigma2.alpha * (ybar.i - mu) ) / (1/var)
    mean <- ( J/sigma2.epsilon * (ybar.i - mu) ) / (1/var)
    #mean <- ( J*sigma2.alpha * (ybar.i - mu) ) / (1/var)
    #mean <- ( J/sigma2.epsilon * (ybar.i - mu) ) / (var)
    rnorm(size, mean, sqrt(var))
}

sample.mu.alpha <- function(size, ybar.., alphabar, I, J, sigma2.epsilon) {
    rnorm(size, ybar.. - alphabar, sqrt(sigma2.epsilon/(I*J)))
}  

sample.mu <- function(size, eta, I, J, sigma2.alpha) {
    stopifnot(length(eta) == I)
    rnorm(size, sum(eta)/I, sqrt(sigma2.alpha/J))
}

sample.eta.i <- function(size, ybar.i, sigma2.alpha) {
    rnorm(size, ybar.i, sqrt(sigma2.alpha))
}

sample.eta.i2 <- function(size, y, ybar.., I, J, sigma2.epsilon, sigma2.alpha) {
    var <- 1 / ( J/sigma2.epsilon + 1/sigma2.alpha )
    mean <- (J*ybar.i/sigma2.epsilon + 1/sigma2.alpha) / var
    #mean <- ( (sum(y - ybar..))/sigma2.epsilon + 1/sigma2.alpha) / ( 1/var )
    rnorm(size, mean, sqrt(var))
}

sample.eta2 <- function(size, Y, ybar.i, ybar.., I, J, sigma2.epsilon, sigma2.alpha) {
    stopifnot(length(ybar.i) == I) 
    sapply(1:length(ybar.i), function(i)
        sample.eta.i(size, Y[, i], ybar.., I, J, sigma2.epsilon, sigma2.alpha))
}

mu.params <- function(ybar.i, J, sigma2.epsilon, sigma2.alpha) {
    V.mu <- 1 / ( sum(1/(sigma2.epsilon/J + sigma2.alpha)) )
    mu.hat <- sum(1/ ( (sigma2.epsilon/J + sigma2.alpha)*ybar.i )) / ( 1/V.mu )
    return(cbind(mu.hat = mu.hat, V.mu = V.mu))
}

theta.params <- function(ybar.i, mu, J, sigma2.epsilon, sigma2.alpha) {
    V <- 1 / (J/sigma2.epsilon + 1/sigma2.alpha)
    theta.j <- ( (J/sigma2.epsilon)*ybar.i + (1/sigma2.alpha)*mu ) / ( 1/V )
    return(list(theta.j = theta.j, V = V))
}

sample.theta <- function(size, theta.i, V) {
    sapply(1:length(theta.i), function(i) rnorm(size, theta.i[i], V))
}

model.oneway.anova.theta <- function(Y, init, sigma2.epsilon, sigma2.alpha, iters) {
    ybar.i <- colMeans(Y)
    I <- ncol(Y)
    Phi <- matrix(rep(NA, (iters+1)*(I+1)), nrow = iters+1, ncol = I+1)
    Phi[1, ] <- c(init["mu"], rep(init["alpha.i"], I))
    for (i in 2:(iters+1)) {
        mu <- mu.params(ybar.i, J, sigma2.epsilon, sigma2.alpha)
        mu.new <- rnorm(1, mean = mu[1,1], sd = sqrt(mu[1,2]))        
        theta <- theta.params(ybar.i, mu.new, J, sigma2.epsilon, sigma2.alpha)       
        theta.new <- sample.theta(1, theta$theta.j, theta$V)
        Phi[i, ] <- c(mu.new, theta.new)
    }
    return(Phi[-1, ])    
}
