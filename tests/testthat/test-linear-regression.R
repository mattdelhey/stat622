set.seed(1)
data <- structure(c(-0.87, -10.74, -3.27, -1.97, 7.5, -7.25, 17.05, 4.96, 
10.4, 11.05, 0.26, 2.51, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 23, 22, 22, 25, 27, 20, 31, 
23, 27, 28, 22, 24, 0, 0, 0, 0, 0, 0, 31, 23, 27, 28, 22, 24), .Dim = c(12L, 
5L), .Dimnames = list(NULL, c("uptake", "intercept", "aerobic", 
         "age", "aerobic.age")))

y <- as.matrix(data[, 1])
X <- data[, -1]

###
### OLS
###
f <- ols.regression(X, y)
f <- lapply(f, round, 2)

# See if fit equals Hoff truth
truth <- list(
    beta.hat = c(-51.29, 13.11, 2.09, -0.32)
  , beta.stderr = c(12.25, 15.76, 0.53, 0.65)
  , sigma2.hat = 8.54)

lapply(1:length(f), function(i) {
    expect_equal(as.numeric(f[i][[1]]), as.numeric(truth[i][[1]]))
})

# See if we agree with lm
lm <- lm(uptake ~ . - 1, data = as.data.frame(data))
lm.b <- coef(summary(lm))[, "Estimate"]
lm.be <- coef(summary(lm))[, "Std. Error"]
lm.s <- summary(lm)$sigma^2
lm.l <- list(lm.b, lm.be, lm.s)
lm.l <- lapply(lm.l, round, 2)

lapply(1:length(f), function(i) {
    expect_equal(as.numeric(f[i][[1]]), as.numeric(lm.l[[i]]))
})

###
### gprior Bayesian
###
truthb <- list(
    e.beta = c(-47.35, 12.10, 1.93, -0.29)
  , beta.stderr = c(14.41, 18.62, 0.62, 0.77))

prior <- c(g = 12, nu.0 = 1, sigma2.0 = 8.54)
fb <- gprior.regression(X, y, 1000, prior)

fbc <- list(fb$e.beta, apply(fb$beta, 2, sd))
fbc <- lapply(fbc, round, 2)

lapply(1:length(fbc), function(i) {
    expect_equal(as.numeric(fbc[[i]]), as.numeric(truthb[i][[1]]))
})
