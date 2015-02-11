set.seed(1)

data <- c(2.3656491, 2.4952035, 1.0837817, 0.7586751, 0.8780483, 1.2765341,
          1.4598699, 0.1801679, -1.0093589, 1.4870201, -0.1193149, 0.2578262)

phi <- model.normal.nonconjugate.imh(
    phi.0 = c(mean(data), var(data))
  , imh.samples = 100000
  , imh.burnin = 50000
  , alpha = 2
  , beta = 2
  , m = 1
  , v = 10
  , y = data
    )

mean(phi$theta)
mean(phi$sigma2)

plot.mcmc.lines(phi$theta, phi$burnin)
plot.mcmc.lines(phi$sigma2, phi$burnin)

# Posterior probability that mu is bigger than 0.5
mean(phi$theta > 0.5)

