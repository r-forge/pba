\name{pbaIterateConfounding}
\alias{pbaIterateConfounding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaIterateConfounding(model, bias.tables, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{bias.tables}{
%%     ~~Describe \code{bias.tables} here~~
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
function (model, bias.tables, iter) 
{
    models.updated <- list()
    model.updated <- model
    outcome <- names(model$model)[1]
    for (i in 1:iter) {
        for (i in names(bias.tables)) {
            if (length(bias.tables[[i]]$confounding) > 0) {
                for (j in names(bias.tables[[i]]$confounding)) {
                  exposure <- i
                  confounding <- bias.tables[[i]]$confounding[[j]]
                  model.updated <- pbaIterateConfoundingSingle(exposure = exposure, 
                    outcome = outcome, model = model.updated, 
                    confounding = confounding, name = j, iter = 1)[[1]]
                }
            }
        }
        models.updated[[i]] <- model.updated
    }
    return(models.updated)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
