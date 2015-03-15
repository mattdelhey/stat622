devtools::load_all("~/stat622")
I <- 5
J <- 2

sigma2.epsilon <- 1
sigma2.alpha <- 1
mu.true <- 5

Y <- sample.yij.unconditional(I, J, mu.true, sigma2.epsilon, sigma2.alpha)

init <- c(mu = 0, alpha.i = 0)
Phi <- model.oneway.anova.theta(Y, init, sigma2.epsilon, sigma2.alpha, iters = 1000)

colMeans(Phi)
# OLD: [1]  0.09459396 41.49843691 41.73595235 42.06857533 41.58875780 40.83936724
# NEW: [1] 5.6463893 1.5769434 2.8259697 1.1534403 0.4769508 1.3968843

Phi <- model.oneway.anova.eta(Y, init, sigma2.epsilon, sigma2.alpha, iters = 1000)
colMeans(Phi)
