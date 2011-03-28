\name{rtrapezoid}
\alias{rtrapezoid}
\title{Random generation from a trapezoidal distribution.}
\usage{rtrapezoid(n, min=0, mode1=0.33, mode2=0.67, max=1)}
\description{Random generation from a trapezoidal distribution.}
\value{A vector numbers randomly generated from the trapezoidal distribution.}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\references{Adapted from Matthew Fox and colleagues' SAS macro SENSITIVITY ANALYSIS MISCLASSIFICATION MACRO version 1.1.
Available: \link{http://sites.google.com/site/biasanalysis/}}
\arguments{\item{n}{Number of observations. If length(n) > 1, the length is taken to be the number required.}
\item{min}{Lower bound of the trapezoid.}
\item{mode1}{Lower mode of the trapezoid.}
\item{mode2}{Upper mode of the trapezoid.}
\item{max}{Upper bound of the trapezoid.}}
