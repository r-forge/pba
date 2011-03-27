\name{pbaCalculatePredictiveValuesInternal}
\alias{pbaCalculatePredictiveValuesInternal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaCalculatePredictiveValuesInternal(exposure, misclassification, model, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exposure}{
%%     ~~Describe \code{exposure} here~~
}
  \item{misclassification}{
%%     ~~Describe \code{misclassification} here~~
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
function (exposure, misclassification, model, iter) 
{
    table.star <- table(model$model[, 1], model$model[, exposure])
    a1.star <- table.star[4]
    a0.star <- table.star[2]
    b1.star <- table.star[3]
    b0.star <- table.star[1]
    a1.stars <- rep(a1.star, iter)
    a0.stars <- rep(a0.star, iter)
    b1.stars <- rep(b1.star, iter)
    b0.stars <- rep(b0.star, iter)
    result <- pbaBackCalculate(a1.star = a1.stars, a0.star = a0.stars, 
        b1.star = b1.stars, b0.star = b0.stars, se.a = misclassification$se.a, 
        sp.a = misclassification$sp.a, se.b = misclassification$se.b, 
        sp.b = misclassification$sp.b)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
