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
        theta.proposal <- rbeta(1, alpha, beta)

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
        sigma2.proposal <- rlnorm(1, m, sqrt(v))

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
      dlnorm(sigma2.proposal, m, sqrt(v), log = TRUE) 
    
    sigma2.state.posterior <- sum(dnorm(y, theta.state, sqrt(sigma2.state), log = TRUE)) +
      dlnorm(sigma2.state, m, sqrt(v), log = TRUE)

    J.state <- 0
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

    J.state <- 0
    J.proposal <- 0

    log.r <- (theta.proposal.posterior + J.state) - (theta.state.posterior + J.proposal)
    return(log.r)
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
