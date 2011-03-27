\name{pbaCorrectMisclassification}
\alias{pbaCorrectMisclassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaCorrectMisclassification(bias.tables, pba.variables, model, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bias.tables}{
%%     ~~Describe \code{bias.tables} here~~
}
  \item{pba.variables}{
%%     ~~Describe \code{pba.variables} here~~
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
function (bias.tables, pba.variables, model, iter) 
{
    for (i in names(pba.variables)) {
        rows <- which(bias.tables[[i]]$misclassification$a1.hat < 
            0 | bias.tables[[i]]$misclassification$a0.hat < 0 | 
            bias.tables[[i]]$misclassification$b1.hat < 0 | bias.tables[[i]]$misclassification$b0.hat < 
            0 | bias.tables[[i]]$misclassification$ppv.a < 0 | 
            bias.tables[[i]]$misclassification$ppv.a > 1 | bias.tables[[i]]$misclassification$npv.a < 
            0 | bias.tables[[i]]$misclassification$npv.a > 1 | 
            bias.tables[[i]]$misclassification$ppv.b < 0 | bias.tables[[i]]$misclassification$ppv.b > 
            1 | bias.tables[[i]]$misclassification$npv.b < 0 | 
            bias.tables[[i]]$misclassification$npv.b > 1 | bias.tables[[i]]$misclassification$se.a < 
            0 | bias.tables[[i]]$misclassification$se.a > 1 | 
            bias.tables[[i]]$misclassification$se.b < 0 | bias.tables[[i]]$misclassification$se.b > 
            1 | bias.tables[[i]]$misclassification$sp.a < 0 | 
            bias.tables[[i]]$misclassification$sp.a > 1 | bias.tables[[i]]$misclassification$sp.b < 
            0 | bias.tables[[i]]$misclassification$sp.b > 1)
        while (length(rows) > 0) {
            bias.tables[[i]]$misclassification[rows, c("se.a", 
                "se.b", "sp.a", "sp.b")] <- pbaSampleMisclassification(misclassification = pba.variables[[i]]$misclassification, 
                iter = length(rows))
            result <- pbaCalculatePredictiveValuesInternal(exposure = i, 
                misclassification = bias.tables[[i]]$misclassification[rows, 
                  ], model = model, iter = length(rows))
            bias.tables[[i]]$misclassification[rows, c("a1.hat", 
                "a0.hat", "b1.hat", "b0.hat", "ppv.a", "npv.a", 
                "ppv.b", "npv.b")] <- result[, c("a1.hat", "a0.hat", 
                "b1.hat", "b0.hat", "ppv.a", "npv.a", "ppv.b", 
                "npv.b")]
            rows <- which(bias.tables[[i]]$misclassification$a1.hat < 
                0 | bias.tables[[i]]$misclassification$a0.hat < 
                0 | bias.tables[[i]]$misclassification$b1.hat < 
                0 | bias.tables[[i]]$misclassification$b0.hat < 
                0 | bias.tables[[i]]$misclassification$ppv.a < 
                0 | bias.tables[[i]]$misclassification$ppv.a > 
                1 | bias.tables[[i]]$misclassification$npv.a < 
                0 | bias.tables[[i]]$misclassification$npv.a > 
                1 | bias.tables[[i]]$misclassification$ppv.b < 
                0 | bias.tables[[i]]$misclassification$ppv.b > 
                1 | bias.tables[[i]]$misclassification$npv.b < 
                0 | bias.tables[[i]]$misclassification$npv.b > 
                1 | bias.tables[[i]]$misclassification$se.a < 
                0 | bias.tables[[i]]$misclassification$se.a > 
                1 | bias.tables[[i]]$misclassification$se.b < 
                0 | bias.tables[[i]]$misclassification$se.b > 
                1 | bias.tables[[i]]$misclassification$sp.a < 
                0 | bias.tables[[i]]$misclassification$sp.a > 
                1 | bias.tables[[i]]$misclassification$sp.b < 
                0 | bias.tables[[i]]$misclassification$sp.b > 
                1)
        }
    }
    return(bias.tables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
