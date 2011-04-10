\name{summary.pba}
\alias{summary.pba}
\title{summary.pba}
\usage{summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
summary.pba(pba, transformation, scale="additive", alpha=0.05, ...)
}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\arguments{\item{pba}{A \code{pba} pbject.}
\item{transformation}{(optional) A character string naming a function to 
apply to the coefficients. If using a function that transforms the 
coefficients to a multiplicative scale, be sure to specify the \code{scale} 
argument accordingly.}
\item{scale}{A character string of either \code{"additive"} or 
\code{"multiplicative"}. If \code{"additive"}, precision will be calculated 
the upper confidience level minus the lower confidence level. If 
\code{"multiplicative"}, the precision will be calculated as the upper 
confidence level divided by the lower confidence interval. Precision is 
calculated after any transformations.}
\item{alpha}{Alpha level used for determining confidence limits.}
\item{...}{(optional) Additional arguments to pass to the 
\code{transformation} function, if specified.}}
\alias{summary.pba}

\alias{summary.pba}

\alias{summary.pba}

\alias{summary.pba}

\alias{summary.pba}

\alias{summary.pba}

