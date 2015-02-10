alpha.1 <- 5
beta.1  <- 2
delta.1 <- 0.3

alpha.2 <- 2
beta.2  <- 8
delta.2 <- 0.7

n.samples <- 100000

mixture <- rbetamixture(
    n.samples
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , delta.1 = delta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta.2 = delta.2)
hist(mixture, prob = TRUE)
lines(density(mixture))

square <- function(x) x^2
square.minus <- function(x, mu) (x - mu)^2

expectation.mu <- model.mixture.beta.importance(
    n.samples
  , h       = identity
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , delta.1 = delta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta.2 = delta.2)

expectation.sigma2 <- model.mixture.beta.importance(
    n.samples
  , h       = square.minus
  , mu      = expectation.mu
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , delta.1 = delta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta.2 = delta.2)

expectation.square <- model.mixture.beta.importance(
    n.samples
  , h       = sqrt
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , delta.1 = delta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta.2 = delta.2)

mean(mixture)
expectation.mu

var(mixture)
expectation.sigma2
expectation.square - expectation.mu^2

mean(mixture >= 0.45 & mixture <= 0.55)

indicator <- function(x, lower, upper) if (x >= lower && x <= upper) x else 0

expectation.prob <- model.mixture.beta.importance(
    n.samples
  , h       = indicator
  , lower   = 0.45
  , upper   = 0.55
  , alpha.1 = alpha.1
  , beta.1  = beta.1
  , delta.1 = delta.1
  , alpha.2 = alpha.2
  , beta.2  = beta.2
  , delta.2 = delta.2)

expectation.prob
