% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABC.R
\name{ABC_on_hyperparameters}
\alias{ABC_on_hyperparameters}
\title{Conduct an ABC on a range of effect size distribution hyperparameters}
\usage{
ABC_on_hyperparameters(
  x,
  phenotypes,
  iters,
  effect_distribution = rbayesB,
  parameter_distributions = list(pi = function(x) rbeta(x, 25, 1), d.f = function(x)
    runif(x, 1, 100), scale = function(x) rbeta(x, 1, 3) * 100),
  h_dist = function(x) rep(0.5, x),
  center = TRUE,
  par = FALSE,
  phased = FALSE,
  save_effects = FALSE
)
}
\arguments{
\item{x}{matrix. Input genotypes, SNPs as rows, columns as individuals. Genotypes formatted as 0,1,2 for the major homozygote, heterozygote, and minor homozygote, respectively.}

\item{phenotypes}{numeric vector. Observed phenotypes, one per individual.}

\item{iters}{numeric. Number of ABC permutations to run.}

\item{effect_distribution}{function, default rbayesB. A function for a distribution of effect sizes. The first argument, n,
should be the number of effects to draw. Other arguments must be named to match the names of the parameter distribution
functions given in the parameter_distributions argument.}

\item{parameter_distributions}{List containing named functions, defualt
list(pi = function(x) rbeta(x, 25, 1), d.f = function(x) runif(x, 1, 100), scale = function(x) rbeta(x, 1, 3)*100).
A list containing functions for each parameter in the effect size distribution. The given functions should be named to
match the parameter for which they produce distributions, and should take a single argument, x, which holds the number of
parameter values to draw.}

\item{h_dist}{function, default function(x) rep(.5, x). A function for the distribution of heritability from which to draw
h^2 values for each iteration. Must take a single argument, x, which holds the number of h values to draw. The default
produces only heritabilities of .5.}

\item{center}{logical, default T. Determines if the phenotypes provided and generated during each iteration
should be centered (have their means set to 0).}

\item{par}{numeric or FALSE, default FALSE. If numeric, the number of cores on which to run the ABC.}
}
\description{
Runs Approximate Bayesian Computation across a range of marker effect size
distribution hyperparameters by comparing the distributions of phenotypes produced to those in the real data.
}
\details{
Please note that no missing phenotypes or missing genotype data is permitted. Missing genotypes should be imputed
(such as with \code{\link{impute_and_phase_beagle}}) before this function is run.
}
