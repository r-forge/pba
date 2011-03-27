\name{pbaBackCalculateConfounding}
\alias{pbaBackCalculateConfounding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaBackCalculateConfounding(a1. = NULL, a0. = NULL, b1. = NULL, b0. = NULL, p1 = NULL, p0 = NULL, rr = NULL, rd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a1.}{
%%     ~~Describe \code{a1.} here~~
}
  \item{a0.}{
%%     ~~Describe \code{a0.} here~~
}
  \item{b1.}{
%%     ~~Describe \code{b1.} here~~
}
  \item{b0.}{
%%     ~~Describe \code{b0.} here~~
}
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{p0}{
%%     ~~Describe \code{p0} here~~
}
  \item{rr}{
%%     ~~Describe \code{rr} here~~
}
  \item{rd}{
%%     ~~Describe \code{rd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (a1. = NULL, a0. = NULL, b1. = NULL, b0. = NULL, p1 = NULL, 
    p0 = NULL, rr = NULL, rd = NULL) 
{
    n11 <- (a1. + b1.) * p1
    n01 <- (a0. + b0.) * p0
    n10 <- (a1. + b1.) - n11
    n00 <- (a0. + b0.) - n10
    n1. <- a1. + b1.
    n0. <- a0. + b0.
    if (!is.null(rr)) {
        a11 <- (rr * n11 * a1.)/(rr * n11 + n1. - n11)
        a01 <- (rr * n01 * a0.)/(rr * n01 + n0. - n01)
        a10 <- a1. - a11
        a00 <- a0. - a01
        b11 <- n11 - a11
        b01 <- n01 - a01
        b10 <- b1. - b11
        b00 <- b0. - b01
    }
    else {
        a11 <- ((rd * n1. * (n1. - n11)) + (a1. * n11))/n1.
        a01 <- ((rd * n0. * (n0. - n01)) + (a0. * n01))/n0.
        a10 <- a1. - a11
        a00 <- a0. - a01
        b11 <- n11 - a11
        b01 <- n01 - a01
        b10 <- b1. - b11
        b00 <- b0. - b01
    }
    result <- data.frame(a11, a01, a10, a00, b11, b01, b10, b00)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
