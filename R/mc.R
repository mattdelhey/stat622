#' @title mc.quantile.ci
#' Monte Carlo quantile-based confidence interval
mc.quantile.ci <- function(posterior.sample, alpha = 0.05) {
    posterior.mean <- mean(posterior.sample)
    posterior.quantiles <- quantile(posterior.sample, c(alpha/2, 1 - alpha/2))
    
    mc.quantile.ci <- list(
        posterior.mean = posterior.mean
      , posterior.qi.lower = posterior.quantiles[1]
      , posterior.qi.upper = posterior.quantiles[2])
    return(mc.quantile.ci)
}
