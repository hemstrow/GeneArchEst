% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs.R
\name{optimum_constant_slide}
\alias{optimum_constant_slide}
\title{Increase the selection optimum by a given percentage of starting variance each generation.}
\usage{
optimum_constant_slide(opt, iv, slide = 0.3, ...)
}
\arguments{
\item{opt}{numeric. Initial optimum phenotype.}

\item{iv}{numeric. Initial variance. The optimum will increase by some portion of this variance.}

\item{slide}{numeric. Proportion of the initial variance by which the optimum phenotype will slide.}
}
\description{
Increase the selection optimum by a given percentage of starting variance each generation.
}
