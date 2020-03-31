#=======distribution functions=========
#' Get random draws from the distribution used for bayesB regressions.
#'
#' Generate any number of random values drawn from the distribution used for
#' bayesB genomic regressions.
#'
#' Under a bayesB model, the effect of each site is drawn from a distribution
#' where var(g) == 0 with prob \emph{pi}, and is otherwise drawn from a scaled
#' t distribution with degrees of freedom \emph{d.f} and scale
#' \emph{scale}.
#'
#' @param n numeric. Number of draws to make from the distribution
#' @param pi numeric. Probability that any one site has zero effect
#' @param d.f numeric. Degrees of freedom for the scaled-inverse chi-squared
#'   distribution
#' @param scale numeric. Scale/shape parameter for the scaled-inverse
#'   chi-squared distribution.
#'
#' @export
rbayesB <- function(n, pi, d.f, scale){
  effects <- rbinom(n, 1, 1 - pi) # these are non-zero
  #effects[effects != 0] <- LaplacesDemon::rinvchisq(sum(effects), d.f, scale) # inverse chi distribution alternative
  effects[effects != 0] <- scale * rt(sum(effects), d.f)
  return(effects)
}

