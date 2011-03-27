\name{pbaPlotSelection}
\alias{pbaPlotSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaPlotSelection(pba = NULL, data = NULL, density = T, parameters = c("s.a1", "s.a0", "s.b1", "s.b0"), title = "Selection", scales = "free", L2 = "selection")
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
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{L2}{
%%     ~~Describe \code{L2} here~~
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
function (pba = NULL, data = NULL, density = T, parameters = c("s.a1", 
    "s.a0", "s.b1", "s.b0"), title = "Selection", scales = "free", 
    L2 = "selection") 
{
    pbaPlotMisclassification(pba = pba, data = data, density = density, 
        parameters = parameters, title = title, scales = scales, 
        L2 = L2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
