#=======distribution functions=========
#' Effect size distributions.
#'
#' Generate any number of random values drawn from the locus effect size
#' distributions.
#'
#' @section BayesB:
#' Under a bayesB model, each locus has probability \emph{pi} of having no effect. Effects are otherwise drawn from a scaled t
#' distribution with degrees of freedom \emph{d.f} and scale \emph{scale}.
#' Alternatively, a fixed number of sites (\emph{sites}) can be given an effect from the same
#' distribution. This is equivalent to the expected number of sites when \emph{pi} = 1
#' - (\emph{sites}/\emph{n}), where \emph{n} is the number of loci.
#'
#' @section BayesA:
#' Under a bayesA model, all loci have effects drawn from a scaled t distribution
#' with degrees of freedom \emph{d.f} and scale \emph{scale}. This is equivalent
#' to a bayesB model where \emph{pi} is equal to 0.
#'
#' @section scaled normal:
#' Effects are drawn from a normal distribution where var(g) == 0.
#' Each locus has probability \emph{pi} of having no effect. Effects are otherwise drawn from a normal
#' distribution with mean 0, standard deviation \emph{sd}, and scale \emph{scale}.
#' Alternatively, a fixed number of sites (\emph{sites}) can be given an effect from the same
#' distribution. This is equivalent to the expected number of sites when \emph{pi} = 1
#' - (\emph{sites}/\emph{n}), where \emph{n} is the number of loci.
#'
#' @section flat:
#' Each locus has probability \emph{pi} of having no effect. Effects are otherwise equal to
#' the parameter scale \emph{scale}.
#' Alternatively, a fixed number of sites (\emph{sites}) can be given an effect from the same
#' distribution. This is equivalent to the expected number of sites when \emph{pi} = 1
#' - (\emph{sites}/\emph{n}), where \emph{n} is the number of loci.
#'
#' @param n numeric. Number of draws to make from the distribution
#' @param pi numeric. Probability that any one site has zero effect
#' @param sites numeric. Number of sites that have an effect.
#' @param d.f numeric. Degrees of freedom for the effect size t distribution.
#' @param sd numeric. Standard deviation for the effect size t distribution.
#' @param scale numeric. Scale/shape parameter for the scaled t distribution.
#' @name effect_size_distributions
#'
#' @aliases rbayesB rbayesB_fixed rzero_inflated_normal rzero_inflated_normal_fixed rflat rflat_fixed
NULL



#' @describeIn effect_size_distributions assign some loci an effect from a scaled t distribution given a probablity of non-zero effects
#' @export
rbayesB <- function(n, pi, d.f, scale){
  eff <- rbinom(n, 1, 1 - pi) # these are non-zero
  #eff[eff != 0] <- LaplacesDemon::rinvchisq(sum(eff), d.f, scale) # inverse chi distribution alternative
  eff[eff != 0] <- scale * rt(sum(eff), d.f)
  return(eff)
}

#' @describeIn effect_size_distributions assign a fixed number of loci effects from a scaled t distribution
#' @export
rbayesB_fixed <- function(n, sites, d.f, scale){
  sites <- floor(sites)
  eff <- numeric(n)
  eff[sample(n, sites, replace = F)] <- scale * rt(sites, d.f)
  return(eff)
}

#' @describeIn  effect_size_distributions assign all loci effects from a scaled t distribution
#' @export
rbayesA <- function(n, d.f, scale){
  eff <- scale*rt(n, d.f)
}

#' @describeIn effect_size_distributions assign some loci an effect from a scaled normal distribution given a probability of non-zero effects
#' @export
rzero_inflated_normal <- function(n, pi, sd, scale){
  eff <- rbinom(n, 1, 1 - pi) #does each site have an effect?
  eff[which(eff == 1)] <- scale * rnorm(sum(eff), 0, sd) #what are the effect sizes?
  return(eff)
}

#' @describeIn effect_size_distributions assign a fixed number of loci effects from a scaled normal distribution
#' @export
rzero_inflated_normal_fixed <- function(n, sites, sd, scale){
  sites <- floor(sites)
  eff <- numeric(n)
  eff[sample(n, sites, replace = F)] <- scale * rnorm(sites, 0, sd)
  return(eff)
}

#' @describeIn effect_size_distributions asign loci a fixed effect given a probability of non-zero effects
#' @export
rflat <- function(n, pi, scale){
  eff <- rbinom(n, 1, 1 - pi) #does each site have an effect?
  eff[which(eff == 1)] <- scale
  return(eff)
}

#' @describeIn effect_size_distributions assign a fixed number of loci a fixed effect
#' @export
rflat_fixed <- function(n, sites, scale){
  sites <- floor(sites)
  eff <- numeric(n)
  eff[sample(n, sites, replace = F)] <- scale
  return(eff)
}


#' Reasonable parameter transformations
#'
#' Produce reasonable linear transformations for any provided parameters from a
#' bank of tested transformation functions. If no tested transforms exist, returns
#' the \code{\link{identity}} function.
#'
#' @param parameters character, default c("pi", "scale"). Parameters for which to fetch reasonable transformation functions.
#'
#' @return A nested, named list containing forward and back transforms for the named parameters.
#'
#' @export
#'
#' @examples
#' # for pi
#' reasonable_transform("pi")
#'
#' # for pi and scale
#' reasonable_transform(c("pi", "scale"))
#'
#' # when there aren't any banked transformations, returns the identity function
#' reasonable_transform(c("pi", "not_a_parameter", "scale"))
reasonable_transform <- function(parameters = c("pi", "scale")){
  #===========bank of often reasonable transforms===
  # forward transforms
  forwards <- list(pi = function(pi) log10(1 - pi),
                   scale = function(scale) log10(scale),
                   sites = function(sites) log10(sites))

  # back transforms
  back <- list(pi = function(pi) 1 - 10^pi,
               scale = function(scale) 10^scale,
               sites = function(sites) 10^sites)



  #===========get the selected transforms============
  # initialize
  st <- vector("list", 2)
  names(st) <- c("forwards", "back")
  for(i in 1:length(st)){
    st[[i]] <- vector("list", length(parameters))
    names(st[[i]]) <- parameters
  }

  # return the banked transform if possible, otherwise return an identity function
  for(i in parameters){
    st[[1]][[i]] <- ifelse(i %in% names(forwards), forwards[[i]], identity)
    st[[2]][[i]] <- ifelse(i %in% names(back), back[[i]], identity)
  }

  return(st)
}
