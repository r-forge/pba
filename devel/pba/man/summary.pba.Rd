\name{summary.pba}
\alias{summary.pba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.pba(pba, transformation = NULL, scale = "additive", alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pba}{
%%     ~~Describe \code{pba} here~~
}
  \item{transformation}{
%%     ~~Describe \code{transformation} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (pba, transformation = NULL, scale = "additive", alpha = 0.05, 
    ...) 
{
    model <- pba$model
    coefficients.star <- pba$coefficients.star
    coefficients.hat <- pba$coefficients.hat
    coefficients.hat.random <- pba$coefficients.hat.random
    summary <- list()
    estimate <- coefficients.star[, "Estimate"]
    summary$star <- data.frame(estimate = coefficients.star[, 
        "Estimate"], confint(model, level = 1 - alpha), p = coefficients.star[, 
        "Pr(>|z|)"])
    names(summary$star) <- c("estimate", paste(100 * alpha/2, 
        "\%", sep = ""), paste(100 * (1 - alpha/2), "\%", sep = ""), 
        "Pr(>|z|)")
    summary$hat <- t(sapply(coefficients.hat, function(x) {
        quantile(x[, "Estimate"], c(0.5, alpha/2, 1 - alpha/2))
    }))
    summary$hat <- data.frame(summary$hat)
    names(summary$hat) <- c("estimate", paste(100 * alpha/2, 
        "\%", sep = ""), paste(100 * (1 - alpha/2), "\%", sep = ""))
    summary$hat$p <- as.vector(sapply(coefficients.hat, function(x) {
        length(x[, "Estimate"][x[, "Estimate"] < 0])/length(x[, 
            "Estimate"])
    }))
    summary$hat$p <- ifelse(summary$hat$p > 0.5, 1 - summary$hat$p, 
        summary$hat$p)
    summary$hat.random <- t(sapply(coefficients.hat.random, function(x) {
        quantile(x, c(0.5, alpha/2, 1 - alpha/2))
    }))
    summary$hat.random <- data.frame(summary$hat.random)
    names(summary$hat.random) <- c("estimate", paste(100 * alpha/2, 
        "\%", sep = ""), paste(100 * (1 - alpha/2), "\%", sep = ""))
    summary$hat.random$p <- as.vector(sapply(coefficients.hat.random, 
        function(x) {
            length(x[x < 0])/length(x)
        }))
    summary$hat.random$p <- ifelse(summary$hat.random$p > 0.5, 
        1 - summary$hat.random$p, summary$hat.random$p)
    if (!is.null(transformation)) {
        summary <- lapply(summary, function(x) {
            x[, -4] <- do.call(transformation, args = list(x[, 
                -4], ...))
            return(x)
        })
    }
    if (scale == "additive") {
        summary <- lapply(summary, function(x) {
            x$precision <- x[, 3] - x[, 2]
            return(x[, c(1, 2, 3, 5, 4)])
        })
    }
    if (scale == "multiplicative") {
        summary <- lapply(summary, function(x) {
            x$precision <- x[, 3]/x[, 2]
            return(x[, c(1, 2, 3, 5, 4)])
        })
    }
    class(summary) <- "summary.pba"
    return(summary)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
