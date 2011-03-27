\name{plot.pba}
\alias{plot.pba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.pba(pba = NULL, data = NULL, density = T, exp = F, adjust = 1, binwidth = NULL, scales = "free", variables = NULL, print = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pba}{
%%     ~~Describe \code{pba} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{density}{
%%     ~~Describe \code{density} here~~
}
  \item{exp}{
%%     ~~Describe \code{exp} here~~
}
  \item{adjust}{
%%     ~~Describe \code{adjust} here~~
}
  \item{binwidth}{
%%     ~~Describe \code{binwidth} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{variables}{
%%     ~~Describe \code{variables} here~~
}
  \item{print}{
%%     ~~Describe \code{print} here~~
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
function (pba = NULL, data = NULL, density = T, exp = F, adjust = 1, 
    binwidth = NULL, scales = "free", variables = NULL, print = T) 
{
    plotEstimates(pba = pba, data = data, density = density, 
        exp = exp, adjust = adjust, binwidth = binwidth, scales = scales, 
        variables = variables, print = print)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
