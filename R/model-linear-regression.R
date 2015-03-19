ols.regression <- function(X, y, add.intercept = FALSE) {
    stopifnot(is.matrix(X), is.matrix(y))

    if (add.intercept) {
        X <- add.intercept(X)
    }
    
    beta.hat <- ols.beta.hat(X, y)
    ssr <- ssr.beta(X, y, beta.hat)
    sigma2.hat <- ssr / (nrow(X) - ncol(X))
    beta.var <- ols.beta.var(X, as.numeric(sigma2.hat))
    beta.stderr <- sqrt(diag(beta.var))

    list(
        beta.hat    = beta.hat
      , beta.stderr = beta.stderr
      , sigma2.hat  = sigma2.hat)
}

ols.beta.hat <- function(X, y) {
    solve(t(X) %*% X) %*% t(X) %*% y
}

ols.beta.var <- function(X, sigma2.hat) {
    solve(t(X) %*% X) * sigma2.hat
}

ssr.beta <- function(X, y, b) {
    t(y) %*% y - 2%*%t(b)%*%t(X)%*%y + t(b)%*%t(X)%*%X%*%b
}

add.intercept <- function(X) {
    cbind(intercept = 1, X)
}

gprior.regression <- function(X, y, size, prior, add.intercept = FALSE) {    
    stopifnot(is.matrix(X), is.matrix(y), is.vector(prior))
    
    if (add.intercept) {
        X <- add.intercept(X)
    }

    H <- gprior.h(prior["g"], X)
    ssr.g <- gprior.ssr(y = y, g = prior["g"], H = H)
    sigma2 <- sample.sigma2(size, length(y), ssr.g, prior["nu.0"], prior["sigma2.0"])

    vb <- gprior.vs(prior["g"], X)
    eb <- gprior.ms(vb, X, as.numeric(y))
    
    beta <- sample.mvn(eb, vb, sigma2)
    e.beta <- expectation.beta(prior["g"], ols.beta.hat(X, y))

    f <- list(beta = beta, sigma2 = sigma2, e.beta = e.beta)
    class(f) <- "reg"
    f
}

predict.regression <- function(beta.hat, Xnew, add.intercept = FALSE) {
    if (add.intercept) {
        Xnew <- add.intercept(Xnew)
    }
    
    Xnew %*% beta.hat
}

plot.regression <- function(y, yhat) {
    mse <- mse(y, yhat) #mad(y,yhat)
    plot(y, yhat, main = sprintf("mse: %.3f", mse))
    mse
}

print.reg <- function(f, ...) {    
    x <- as.matrix(cbind(
        f$e.beta
      , apply(f$beta, 2, sd)
      , t(apply(f$beta, 2, quantile, c(0.025, 0.975)))
        ))
    colnames(x) <- c("expectation", "std err", "0.25%", "0.975%")
    
    zero <- abs(sign(x[,3]) + sign(x[,4])) == 2
    x <- cbind(x, zero)
    printCoefmat(x)
}

mse2 <- function(y, yhat) {
    1/length(y) * sum( (y - yhat)^2 )
}

mse <- function(y, yhat) {
    mean( (y - yhat)^2 )
}

mad <- function(y, yhat) {
    mean( abs(y - yhat))
}


expectation.beta <- function(g, ols.beta.hat) {
    (g/(g+1)) * ols.beta.hat
}
    
sample.mvn <- function(eb, vb, sigma2) {
    p <- length(eb); size <- length(sigma2)
    e <- matrix(rnorm(size*p, 0, sqrt(sigma2)), size, p)
    t(t(e %*% chol(vb)) + c(eb))
}

sample.sigma2 <- function(size, n, ssr.g, nu.0, sigma2.0) {
    1/rgamma(size, (nu.0 + n)/2, (nu.0*sigma2.0 + ssr.g)/2)
}

gprior.h <- function(g, X) {
    (g/(g+1)) * X %*% solve(t(X) %*% X) %*% t(X)
}

gprior.ssr <- function(X = NULL, y, g, H = NULL) {
    if (is.null(X)) {
        return(t(y) %*% ( diag(1, nrow = length(y)) - H ) %*% y)
    }
    if (is.null(H)) {
        return(t(y) %*% ( diag(1, nrow(X)) - (g/(g+1))* X%*%solve(t(X) %*% X) %*% t(X) )*y)
    }    
}
    
gprior.ssr.mv <- function(y, m, V) {    
    t(y) %*% y - t(m) %*% solve(V) %*% m
}

gprior.vs <- function(g, X) {
    (g/(g+1)) * solve(t(X) %*% X)
}

gprior.ms <- function(vs, X, y) {
    vs %*% t(X) %*% y
}

gprior.v <- function(g, sigma2, X) {
    (g/(g+1)) * sigma2 * solve(t(X) %*% X)
}

gprior.m <- function(g, X, y) {
    (g/(g+1)) * solve(t(X) %*% X) %*% t(X) * y
}


#' Split data into arbitrary partitions
#' @param n number of indicies
#' @param p vector that sums to one
#' @return List of split indicies
#' @family utils
#' @export
split.data <- function(n, p) {
    if (sum(p) != 1) stop("Split proportions must sum to one.")
    if (length(n) != 1) stop("n is the number of indicies")

    k <- length(p)
    s <- sapply(2:k, function(j) floor(p[j] * n))
    s <- c(n - sum(s), s)
    stopifnot(sum(s) == n)
    
    ind <- vector("list", k)
    ind[[1]] <- sample(n, size = s[1], replace = FALSE)
    for (j in 2:k) {
        available <- setdiff(1:n, unique(unlist(ind[1:j])))
        ind[[j]] <- sample(available, size = s[j], replace = FALSE)
    }
    stopifnot(all(sort(unlist(ind)) == 1:n))

    names(ind) <- names(p)
    ind
}
