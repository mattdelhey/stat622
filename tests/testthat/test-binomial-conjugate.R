a<-1  ; b<-1   #prior
n<-10 ; y<-2   #data
(analytic <- qbeta(c(.025,.975), a+y,b+n-y))

posterior <- model.binomial.conjugate(n.obs = n, y.sum = y, alpha = a, beta = b)
theta.mc <- model.binomial.conjugate.mc(posterior, 10000)$theta.posterior.sample
(mc <- mc.quantile.ci(theta.mc))

expect_equal(analytic[1], as.numeric(mc[[2]]), tolerance = 0.01)
expect_equal(analytic[2], as.numeric(mc[[3]]), tolerance = 0.01)
