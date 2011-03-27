\name{pba}
\alias{pba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pba(model, pba.variables, iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
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
function (model, pba.variables, iter = 1000) 
{
    time.start <- Sys.time()
    call <- match.call()
    model.summary <- summary(model)
    model.summaries <- list()
    bias.lists <- list()
    for (i in 1:iter) {
        if (i\%\%10 == 0) 
            print(i)
        bias.tables <- pbaBiasTables(pba.variables, iter = 1)
        bias.tables <- pbaCalculatePredictiveValues(bias.tables = bias.tables, 
            model = model, iter = 1)
        model.updated <- model
        bias.tables <- pbaCorrectMisclassification(bias.tables = bias.tables, 
            pba.variables = pba.variables, model = model, iter = 1)
        model.updated <- pbaIterateMisclassification(model = model.updated, 
            bias.tables = bias.tables, iter = 1)[[1]]
        model.updated <- pbaIterateSelection(model = model.updated, 
            bias.tables = bias.tables, iter = 1)[[1]]
        model.updated <- pbaIterateConfounding(model = model.updated, 
            bias.tables = bias.tables, iter = 1)[[1]]
        model.updated <- update(model.updated, data = model.updated$data, 
            formula = model.updated$formula)
        model.summaries[[i]] <- summary(model.updated)["coefficients"]
        for (j in names(pba.variables)) {
            bias.lists[[j]][[i]] <- bias.tables[[j]]
        }
    }
    bias.tables <- pbaBiasListsTables(bias.lists = bias.lists)
    coefficients.hat <- list()
    for (i in rownames(model.summaries[[1]]$coefficients)) {
        coefficients.hat[[i]] <- t(sapply(model.summaries, function(x) {
            x$coefficients[i, ]
        }))
    }
    coefficients.hat.random <- list()
    for (i in names(coefficients.hat)) {
        coefficients.hat.random[[i]] <- coefficients.hat[[i]][, 
            "Estimate"] + rnorm(nrow(coefficients.hat[[i]]), 
            0, 1) * rep(model.summary$coefficients[, "Std. Error"], 
            nrow(coefficients.hat[[i]]))
    }
    coefficients.star <- summary(model)$coefficients
    time.end <- Sys.time()
    time.elapsed <- time.end - time.start
    results <- list(call = call, bias.tables = bias.tables, coefficients.hat = coefficients.hat, 
        coefficients.hat.random = coefficients.hat.random, coefficients.star = coefficients.star, 
        bias.lists = bias.lists, model = model, iter = iter, 
        pba.variables = pba.variables, model.summaries = model.summaries, 
        time.elapsed = time.elapsed)
    class(results) <- "pba"
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
