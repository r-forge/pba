\name{pbaSampleConfounding}
\alias{pbaSampleConfounding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaSampleConfounding(confounding, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{confounding}{
%%     ~~Describe \code{confounding} here~~
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
function (confounding, iter) 
{
    results <- list()
    for (i in confounding) {
        if (!is.null(i$p1.distr)) {
            i <- pbaAddIter(i, iter)
            p1 <- do.call(i$p1.distr$distr, i$p1.distr$args)
            p0 <- do.call(i$p1.distr$distr, i$p1.distr$args)
            p1 <- pbaCorrectConfounding(p = p1, distr = i$p1.distr$distr, 
                args = i$p1.distr$args)
            p0 <- pbaCorrectConfounding(p = p0, distr = i$p0.distr$distr, 
                args = i$p0.distr$args)
            if (!is.null(i$rr.distr)) {
                rr <- do.call(i$rr.distr$distr, i$rr.distr$args)
                rd <- rep(NA, iter)
            }
            else {
                rd <- do.call(i$rd.distr$distr, i$rd.distr$args)
                rr <- rep(NA, iter)
            }
        }
        results[[i$name]] <- data.frame(p1, p0, rr, rd)
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
