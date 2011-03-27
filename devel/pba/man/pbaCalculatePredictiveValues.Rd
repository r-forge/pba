\name{pbaCalculatePredictiveValues}
\alias{pbaCalculatePredictiveValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaCalculatePredictiveValues(bias.tables, model, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bias.tables}{
%%     ~~Describe \code{bias.tables} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (bias.tables, model, iter) 
{
    for (i in names(bias.tables)) {
        if (length(bias.tables[[i]]$misclassification) > 0) {
            result <- pbaCalculatePredictiveValuesInternal(exposure = i, 
                misclassification = bias.tables[[i]]$misclassification, 
                model = model, iter = iter)
            bias.tables[[i]]$misclassification <- cbind(bias.tables[[i]]$misclassification, 
                result[, c("a1.hat", "a0.hat", "b1.hat", "b0.hat", 
                  "ppv.a", "npv.a", "ppv.b", "npv.b")])
        }
    }
    return(bias.tables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
