\name{pbaBackCalculate}
\alias{pbaBackCalculate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaBackCalculate(a1.star, a0.star, b1.star, b0.star, se.a, sp.a, se.b, sp.b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a1.star}{
%%     ~~Describe \code{a1.star} here~~
}
  \item{a0.star}{
%%     ~~Describe \code{a0.star} here~~
}
  \item{b1.star}{
%%     ~~Describe \code{b1.star} here~~
}
  \item{b0.star}{
%%     ~~Describe \code{b0.star} here~~
}
  \item{se.a}{
%%     ~~Describe \code{se.a} here~~
}
  \item{sp.a}{
%%     ~~Describe \code{sp.a} here~~
}
  \item{se.b}{
%%     ~~Describe \code{se.b} here~~
}
  \item{sp.b}{
%%     ~~Describe \code{sp.b} here~~
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
function (a1.star, a0.star, b1.star, b0.star, se.a, sp.a, se.b, 
    sp.b) 
{
    b.fp <- 1 - sp.b
    b.fn <- 1 - se.b
    a.fp <- 1 - sp.a
    m1 <- a1.star + a0.star
    m0 <- b1.star + b0.star
    b1.hat <- (sp.b * b1.star - b.fp * b0.star)/(se.b * sp.b - 
        b.fn * b.fp)
    b0.hat <- m0 - b1.hat
    a1.hat <- (a1.star - a.fp * m1)/(se.a + sp.a - 1)
    a0.hat <- m1 - a1.hat
    ppv.a <- (se.a * a1.hat)/((se.a * a1.hat) + ((1 - sp.a) * 
        a0.hat))
    npv.a <- (sp.a * a0.hat)/((sp.a * a0.hat) + ((1 - se.a) * 
        a1.hat))
    ppv.b <- (se.b * b1.hat)/((se.b * b1.hat) + ((1 - sp.b) * 
        b0.hat))
    npv.b <- (sp.b * b0.hat)/((sp.b * b0.hat) + ((1 - se.b) * 
        b1.hat))
    result <- data.frame(a1.hat, a0.hat, b1.hat, b0.hat, se.a, 
        sp.a, se.b, sp.b, ppv.a, npv.a, ppv.b, npv.b)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
