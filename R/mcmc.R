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

#' @title effective.sample.size
#' Calculate effective sample size
effective.sample.size <- function(mcmc.samples, n.lags = 100) {
    n <- length(mcmc.samples)
    ro <- acf(mcmc.samples, lag.max = max(n, n.lags))$acf
    ess <- n / (1 + 2 * sum(ro))
    f <- spectrum(mcmc.samples, method = "ar")    
}

#' @title mc.quantile.ci
#' Monte Carlo quantile-based confidence interval
mc.quantile.ci <- function(posterior.sample, alpha = 0.05) {
    posterior.mean <- mean(posterior.sample)
    posterior.quantiles <- quantile(posterior.sample, c(alpha/2, 50, 1 - alpha/2))
    
    mc.quantile.ci <- list(
        posterior.qi.lower = posterior.quantiles[1]
      , posterior.mean     = posterior.mean
      , posterior.median   = posterior.quantiles[2]
      , posterior.qi.upper = posterior.quantiles[3])
    return(mc.quantile.ci)
}


plot.mcmc.trace <- function(theta, burnin, strip1 = TRUE,
                            xlab = "index", ylab = "theta",
                            ...) {
    if (!require(ggplot2)) stop("ggplot2 required")    

    if (strip1) {
        theta <- theta[-1]
        burnin <- burnin[-1]
    }
    
    qplot(x = 1:length(theta), y = theta, color = burnin, geom = "line", ...) +
      theme.tufte() + xlab(xlab) + ylab(ylab)
}

plot.mcmc.trace.params <- function(theta1, theta2, burnin, strip1 = TRUE,
                                   xlab = "theta1", ylab = "theta2") {
    library(ggplot2)    
    if (length(theta1) != length(theta2))
        stop("theta1 and theta2 must be same length")

    if (strip1) {
        theta1 <- theta1[-1]
        theta2 <- theta2[-1]
        burnin <- burnin[-1]
    }

    qplot(x = theta1, y = theta2, color = burnin, geom = "line") +
      theme.tufte() + xlab(xlab) + ylab(ylab)
}

plot.mcmc.marginal <- function(theta, burnin, strip1 = TRUE,
                               xlab = "theta", ylab = "probability") {
    hist(theta[burnin == "sample"], main = "", xlab = xlab, prob = TRUE)
    lines(density(theta[burnin == "sample"]), col = "blue")
}

plot.mcmc.marginal2 <- function(theta, burnin, strip1 = TRUE,
                               xlab = "theta", ylab = "probability") {
    library(ggplot2)
    # Discard burnin
    theta <- theta[burnin == "sample"]
    breaks <- hist(theta, prob = TRUE)$breaks
    ggplot(data = data.frame(), aes(x = theta)) +
      geom_histogram(aes(y = ..density..), fill = "white", color = "black", breaks = breaks) +
                     #binwidth = (range(theta)[2]-range(theta)[1]) / 100) +
      geom_density(color = "blue") + xlab(xlab) + ylab(ylab) + theme.tufte()
}

    
theme.tufte <- function(base_size = 11, base_family = "serif", ticks = TRUE) {
    # Taken from ggthemes/tufte.R
    ret <- theme_bw(base_family=base_family, base_size=base_size) +
        theme(
            legend.background = element_blank(),
            legend.key        = element_blank(),
            legend.title      = element_blank(),
            legend.position   = "bottom",
            panel.background  = element_blank(),
            panel.border      = element_blank(),
            strip.background  = element_blank(),
            plot.background   = element_blank(),
            axis.line         = element_blank(),
            panel.grid        = element_blank())
    if (!ticks) {
        ret <- ret + theme(axis.ticks = element_blank())
    }
    ret
}

