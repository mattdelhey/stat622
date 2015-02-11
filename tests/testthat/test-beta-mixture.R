set.seed(1)

alpha.1 <- 5
beta.1  <- 2
alpha.2 <- 2
beta.2  <- 8
delta   <- 0.3

n.samples <- 100000

mixture <- rbetamixture(n.samples, alpha.1, beta.1, alpha.2, beta.2, delta)
hist(mixture, prob = TRUE)
lines(density(mixture))

square <- function(x) x^2
square.minus <- function(x, mu) (x - mu)^2
indicator <- function(x, lower, upper) ifelse(x >= lower & x <= upper, 1, 0)

expectation.mu <- model.mixture.beta.importance(
    n.samples
  , h       = identity
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta   = delta)
expect_equal(mean(mixture), expectation.mu, tolerance = 0.05)

expectation.sigma2 <- model.mixture.beta.importance(
    n.samples
  , h       = square.minus
  , mu      = expectation.mu
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta   = delta)
expect_equal(var(mixture), expectation.sigma2, tolerance = 0.05)

expectation.square <- model.mixture.beta.importance(
    n.samples
  , h       = square
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta   = delta)
expect_equal(var(mixture), expectation.square - expectation.mu^2, tolerance = 0.05)

expectation.prob <- model.mixture.beta.importance(
    n.samples
  , h       = indicator
  , lower   = 0.45
  , upper   = 0.55
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta   = delta)
expect_equal(expectation.prob, mean(mixture >= 0.45 & mixture <= 0.55), tolerance = 0.05)
