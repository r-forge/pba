\name{pbaBiasTables}
\alias{pbaBiasTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaBiasTables(pba.variables, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pba.variables}{
%%     ~~Describe \code{pba.variables} here~~
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
function (pba.variables, iter) 
{
    results <- list()
    for (i in pba.variables) {
        if (!is.null(i$misclassification$se.a.distr)) {
            misclassification <- pbaSampleMisclassification(misclassification = i$misclassification, 
                iter = iter)
        }
        else {
            misclassification = list()
        }
        if (!is.null(i$selection$s.a1.distr)) {
            selection <- pbaSampleSelection(selection = i$selection, 
                iter = iter)
        }
        else {
            selection <- list()
        }
        if (!is.null(i$confounding[[1]]$p1.distr)) {
            confounding <- pbaSampleConfounding(confounding = i$confounding, 
                iter = iter)
        }
        else {
            confounding <- list()
        }
        results[[i$variable]] <- list(misclassification = misclassification, 
            selection = selection, confounding = confounding)
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
