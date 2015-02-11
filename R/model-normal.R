#' @title model.normal.conjugate
#' Posterior parameters for the conjugate normal model
model.normal.conjuage <- function(x, mu.0, sigma2.0, kappa.0, nu.0) {
    data <- model.normal.summary(x)
    
    # Posterior parameters
    nu.n <- nu.0 + data$n.obs
    kappa.n <- kappa.0 + data$n.obs
    mu.n <- (kappa.0*mu.0 + data$n.obs * data$y.bar) / kappa.n
    sigma2.n <- (1/nu.n) * (nu.0 * sigma2.0 + (data$n.obs - 1) *
                              data$s2 + (kappa.0*data$n.obs / kappa.n) * (data$y.bar - mu.0)^2)
    
    # Calculate SS.n two ways, check for equivalence
    ss.n.2 <- nu.n * sigma2.n
    ss.n <- (nu.0*sigma2.0) + ss + (n.obs*kappa.0)/kappa.n * (y.bar - mu.0)^2
    stopifnot(all.equal(ss.n, ss.n.2))
    
    posterior.parameters <- list(
        nu.n     = nu.n
      , kappa.n  = kappa.n
      , mu.n     = mu.n
      , sigma2.n = sigma2.n
      , ss.n     = ss.n)
    return(posterior.parameters)
}

#' @title model.normal.conjugate.mc
#' Monte Carlo samples for the conjugate normal model
model.normal.conjugate.mc <- function(posterior.parameters, mc.samples) {
    sigma2.posterior.sample <- 1 / rgamma(
        mc.samples
      , shape = posterior.parameters$nu.n / 2
      , rate  = posterior.parameters$sigma2.n * posterior.parameters$nu.n / 2)
    
    theta.posterior.sample <- rnorm(
        mc.samples
      , mean = posterior.parameters$mu.n
      , sd   = sqrt(sigma2.posterior.sample / posterior.parameters$kappa.n))
    
    predictive.posterior.sample <- rnorm(
        mc.samples
      , mean = theta.posterior.sample
      , sd   = sqrt(sigma2.posterior.sample))

    posterior.samples <- list(
        sigma2.posterior.sample     = sigma2.posterior.sample
      , theta.posterior.sample      = theta.posterior.sample
      , predictive.posterior.sample = predictive.posterior.sample)
    return(posterior.samples)
}

#' @title model.normal.semiconjugate.gibbs
#' Gibbs samples for the semiconjugate normal model
model.normal.semiconjugate.gibbs <- function(phi.0, gibbs.samples, gibbs.burnin = 0,
                                             mu.0, tau2.0, nu.0, sigma2.0,
                                             y.bar, n.obs, s2) {
    if (gibbs.samples %% 1 != 0 | gibbs.burnin %% 1 != 0)
        stop("number of samples must be an integer")

    gibbs.iters <- gibbs.samples + gibbs.burnin + 1
    phi <- construct.phi(phi0, gibbs.iters, gibbs.burnin, vars = c("theta", "sigma2"))

    # Account for phi.0 by starting at 2
    for (i in 2:gibbs.iters) {
        # Use new sigma2 to calculate theta posterior parameters
        theta.conditional <- model.normal.semiconjugate.theta(
            sigma2 = phi$sigma2[i-1]
          , mu.0 = mu.0
          , tau2.0 = tau2.0
          , n.obs = n.obs
          , y.bar = y.bar)

        # Generate new theta
        phi$theta[i] <- rnorm(1, theta.conditional$mu.n,
                              sqrt(theta.conditional$tau2.n))

        # Use new theta to calculate sigma2 posterior parameters
        sigma2.conditional <- model.normal.semiconjugate.sigma2(
            theta = phi$theta[i]
          , nu.0 = nu.0
          , sigma2.0 = sigma2.0
          , n.obs = n.obs
          , y.bar = y.bar
          , s2 = s2)
        
        # Generate new sigma2
        phi$sigma2[i] <- 1 / rgamma(1, sigma2.conditional$nu.n / 2,
                                sigma2.conditional$nu.n * sigma2.conditional$sigma2.n / 2)
    }
    return(phi)
}

#' @title model.normal.semiconjugate.theta
#' Full conditional distribution for theta (sigma2 known)
model.normal.semiconjugate.theta <- function(sigma2, mu.0, tau2.0, n.obs, y.bar) {
    mu.n <- ( mu.0/tau2.0 + n.obs*y.bar/sigma2 ) / ( 1/tau2.0 + n.obs/sigma2 )
    tau2.n <- 1 / ( 1/tau2.0 + n.obs/sigma2 )
    posterior.parameters.theta <- list(
        mu.n = mu.n
      , tau2.n = tau2.n
        )
    return(posterior.parameters.theta)
}

#' @title model.normal.semiconjugate.sigma2
#' Full conditional distribution for sigma2 (theta known)
model.normal.semiconjugate.sigma2 <- function(theta, nu.0, sigma2.0, n.obs, y.bar, s2) {
    nu.n <- nu.0 + n.obs
    # Use trick to calculate sample variance from summary statistics
    sigma2.n <- (1/nu.n) * (nu.0*sigma2.0 + (n.obs-1)*s2 + n.obs*(y.bar - theta)^2)
    posterior.parameters.sigma2 <- list(
        nu.n = nu.n
      , sigma2.n = sigma2.n
        )
    return(posterior.parameters.sigma2)
}

#' @title model.normal.nonconjugate.imh
#' Independent Metropolis-Hastings algorithm for nonconjugate normal model
model.normal.nonconjugate.imh <- function(phi.0, imh.samples, imh.burnin,
                                          alpha, beta, m, v,
                                          y) {
    if (imh.samples %% 1 != 0) stop("number of samples must be an integer")

    accepts.theta <- 0
    accepts.sigma2 <- 0    
    imh.iters <- imh.samples + imh.burnin + 1
    phi <- construct.phi(phi.0, imh.iters, imh.burnin, vars = c("theta", "sigma2"))
    
    for (i in 2:imh.iters) {
        # ~~~ Update theta
        # Sample from (symmetric??) proposal distribtuion (independent of phi)
        theta.proposal <- runif(1, 0, 1)

        # Compute (log) acceptance ratio (ratio of likilihood * prior)       
        log.r <- model.normal.nonconjugate.theta.r(y = y, alpha = alpha, beta = beta
                                                 , theta.proposal = theta.proposal
                                                 , theta.state = phi$theta[i-1]
                                                 , sigma2.state = phi$sigma2[i-1])
        
        # Make descion to keep theta.proposal with probability min(r,1) 
        accept.result <- accept.ratio(log.r, accepts.theta, theta.proposal, phi$theta[i-1])
        phi$theta[i]  <- accept.result$parameter
        accepts.theta <- accept.result$accepts

        # ~~~ Update sigma
        # Sample from (symmetric??) proposal distribution (independent of phi)
        sigma2.proposal <- rlnorm(1, m, sqrt(v)/4)

        log.r <- model.normal.nonconjugate.sigma2.r(y = y, m = m, v = v
                                                  , sigma2.proposal = sigma2.proposal
                                                  , sigma2.state = phi$sigma2[i-1]
                                                  , theta.state = phi$theta[i])
        
        # Make descion to keep theta.proposal with probability min(r,1) 
        accept.result <- accept.ratio(log.r, accepts.sigma2, sigma2.proposal, phi$sigma2[i-1])
        phi$sigma2[i]  <- accept.result$parameter
        accepts.sigma2 <- accept.result$accepts
    }
    
    imh <- list(
        accepts.ratio.theta  = accepts.theta  / imh.iters
      , accepts.ratio.sigma2 = accepts.sigma2 / imh.iters
      , eff.theta  = as.numeric(coda::effectiveSize(phi$theta))
      , eff.sigma2 = as.numeric(coda::effectiveSize(phi$sigma2))
      , phi = phi)
    return(imh)
}

model.normal.nonconjugate.sigma2.r <- function(y, sigma2.proposal, sigma2.state, theta.state, m, v) {
    sigma2.proposal.posterior <- sum(dnorm(y, theta.state, sqrt(sigma2.proposal), log = TRUE)) +
      dlnorm(sigma2.proposal, m, sqrt(v)/4, log = TRUE)
    
    sigma2.state.posterior <- sum(dnorm(y, theta.state, sqrt(sigma2.state), log = TRUE)) +
      dlnorm(sigma2.state, m, sqrt(v)/4, log = TRUE)

    #J.state <- dlnorm(sigma2.state, m, v, log = TRUE)
    J.state <- 0
    #J.proposal <- dlnorm(sigma2.proposal, m, v, log = TRUE)
    J.proposal <- 0

    log.r <- (sigma2.proposal.posterior + J.state) - (sigma2.state.posterior + J.proposal)
    return(log.r)    
}

#' @title model.normal.nonconjugate.theta.r
model.normal.nonconjugate.theta.r <- function(y, theta.proposal, theta.state, sigma2.state, alpha, beta) {
    theta.proposal.posterior <- sum(dnorm(y, theta.proposal, sqrt(sigma2.state), log = TRUE)) +
      dbeta(theta.proposal, alpha, beta, log = TRUE)
    
    theta.state.posterior <- sum(dnorm(y, theta.state, sqrt(sigma2.state), log = TRUE)) +
      dbeta(theta.state, alpha, beta, log = TRUE)

    #J.state <- dunif(theta.state, 0, 1, log = TRUE)
    #J.state <- dbeta(theta.state, 1, 1, log = TRUE)
    J.state <- 0
    #J.proposal <- dunif(theta.proposal, 0, 1, log = TRUE)
    #J.proposal <- dbeta(theta.proposal, 1, 1, log = TRUE)
    J.proposal <- 0

    log.r <- (theta.proposal.posterior + J.state) - (theta.state.posterior + J.proposal)
    return(log.r)
}

#' @title model.normal.summary
#' Summary statistics for normal data
model.normal.summary <- function(x) {
    y.bar <- mean(x)
    s2 <- var(x)
    n.obs <- length(x)
    ss <- (n.obs - 1) * s2

    data <- list(
        y.bar = y.bar
      , s2    = s2
      , n.obs = n.obs
      , ss    = ss)
    return(data)
}

accept.ratio <- function(r, accepts, theta.proposal, theta.state, log = TRUE) {
    # Uniform for accept/reject
    u <- if (log == TRUE) log(runif(1)) else runif(1)
    if (u < r) {
        accepts <- accepts + 1
        res <- theta.proposal
    } else {
        res <- theta.state  
    }
    return(list(parameter = res, accepts = accepts))
}

construct.phi <- function(phi.0, imh.iters, imh.burnin, vars) {
    stopifnot(length(phi.0) == length(vars), is.vector(vars), is.vector(phi.0))
    
    # Subtract one row because we are going to add phi.0
    # Aadd one coloumn because we are going to add a burnin status
    phi.empty <- as.data.frame(matrix(NA, nrow = imh.iters-1, ncol = length(vars)+1))
    names(phi.empty) <- c(vars, "burnin")    
    phi <- rbind(phi.0, phi.empty)

    # Construct burnin status
    phi$burnin[1] <- "init"
    phi$burnin[2:(imh.burnin+1)] <- "burnin"
    phi$burnin[(imh.burnin+1):nrow(phi)] <- "sample"

    return(phi)
}

