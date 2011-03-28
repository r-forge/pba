\name{pbaDistr}
\alias{pbaDistr}
\title{Function to define distributions for bias parameters.}
\usage{pbaDistr(distr, args)}
\description{Function to define distributions for bias parameters.}
\value{An object of class \code{pba.distr}.}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\arguments{\item{distr}{A character vector of length one naming a random generation 
ditribution.}
\item{args}{Additional arguments to be passed to the distr function. The 
argument \code{n} does not need to be specified, as it will be autumatically 
added when appropriate in subsequent functions.}}
