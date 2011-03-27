\name{pbaSampleSelection}
\alias{pbaSampleSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaSampleSelection(selection, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{
%%     ~~Describe \code{selection} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
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
function (selection, iter) 
{
    results <- c()
    if (!is.null(selection$s.a1.distr)) {
        selection <- pbaAddIter(selection, iter)
        s.a1 <- do.call(selection$s.a1.distr$distr, selection$s.a1.distr$args)
        s.a0 <- do.call(selection$s.a0.distr$distr, selection$s.a0.distr$args)
        s.b1 <- do.call(selection$s.b1.distr$distr, selection$s.b1.distr$args)
        s.b0 <- do.call(selection$s.b0.distr$distr, selection$s.b0.distr$args)
        s.a1 <- pbaCorrectConfounding(p = s.a1, distr = selection$s.a1.distr$distr, 
            args = selection$s.a1.distr$args)
        s.a0 <- pbaCorrectConfounding(p = s.a0, distr = selection$s.a0.distr$distr, 
            args = selection$s.a0.distr$args)
        s.b1 <- pbaCorrectConfounding(p = s.b1, distr = selection$s.b1.distr$distr, 
            args = selection$s.b1.distr$args)
        s.b0 <- pbaCorrectConfounding(p = s.b0, distr = selection$s.b0.distr$distr, 
            args = selection$s.b0.distr$args)
    }
    result <- data.frame(s.a1, s.a0, s.b1, s.b0)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
