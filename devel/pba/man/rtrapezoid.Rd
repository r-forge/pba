\name{rtrapezoid}
\alias{rtrapezoid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rtrapezoid(n, min = 0, mode1 = 0.33, mode2 = 0.67, max = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{min}{
%%     ~~Describe \code{min} here~~
}
  \item{mode1}{
%%     ~~Describe \code{mode1} here~~
}
  \item{mode2}{
%%     ~~Describe \code{mode2} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
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
function (n, min = 0, mode1 = 0.33, mode2 = 0.67, max = 1) 
{
    if (length(n) > 1) 
        n <- length(n)
    if (n < 1 | is.na(n)) 
        stop(paste("invalid argument: n =", n))
    n <- floor(n)
    if (any(is.na(c(min, mode1, mode2, max)))) 
        return(rep(NaN, times = n))
    if (min > max | min > mode1 | min > mode2 | mode1 > max | 
        mode2 > max | mode1 > mode2) 
        return(rep(NaN, times = n))
    if (any(is.infinite(c(min, mode1, mode2, max)))) 
        return(rep(NaN, times = n))
    p <- runif(n)
    r <- (p * (max + mode2 - min - mode1) + (min + mode1))/2
    rs.lt.mode1 <- which(r < mode1)
    r[rs.lt.mode1] <- min + sqrt((mode1 - min) * (2 * r[rs.lt.mode1] - 
        min - mode1))
    rs.gt.mode2 <- which(r > mode2)
    r[rs.gt.mode2] <- max - sqrt(2 * (max - mode2) * (r[rs.gt.mode2] - 
        mode2))
    return(r)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
