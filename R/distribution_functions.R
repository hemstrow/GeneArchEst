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

#' @export
zero_inflated_normal <- function(n){
  eff <- rbinom(n, 1, prob.effect) #does each site have an effect?
  eff[which(eff == 1)] <- rnorm(sum(eff), effect.mean, effect.sd) #what are the effect sizes?
  # (x - min(x))/(max(x) - min(x))
  return(eff)
}

#' @export
fixed_number_normal <- function(n){
  eff <- numeric(n)
  eff[sample(n, n.eff, replace = F)] <- rnorm(n.eff, effect.mean, effect.sd)
  return(eff)
}

#' @export
fixed_number_scaled_t <- function(n, sites, pi, d.f, scale){
  effects <- numeric(n)
  effects[sample(n, sites)] <- 1
  effects[effects == 1] <- scale * rt(sum(effects), d.f)
  return(effects)
}

#' @export
fixed_number_flat <- function(n, sites, scale){
  effects <- numeric(n)
  # effects[sample(n, sites)] <- sample(rep(c(scale, -1 * scale), sites), sites, replace = F) # half are positive effects, half are negative, randomly positioned
  effects[sample(n, sites)] <- scale
  return(effects)
}
