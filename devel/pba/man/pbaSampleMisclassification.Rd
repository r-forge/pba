\name{pbaSampleMisclassification}
\alias{pbaSampleMisclassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaSampleMisclassification(misclassification, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{misclassification}{
%%     ~~Describe \code{misclassification} here~~
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
function (misclassification, iter) 
{
    iter.equaled.one <- if (iter == 1) 
        T
    else F
    iter <- if (iter == 1) 
        2
    else iter
    results <- list()
    misclassification <- pbaAddIter(x = misclassification, iter = iter)
    se.as.uncor <- do.call(misclassification$se.a.distr$distr, 
        misclassification$se.a.distr$args)
    sp.as.uncor <- do.call(misclassification$sp.a.distr$distr, 
        misclassification$sp.a.distr$args)
    se.bs.uncor <- do.call(misclassification$se.b.distr$distr, 
        misclassification$se.b.distr$args)
    sp.bs.uncor <- do.call(misclassification$sp.b.distr$distr, 
        misclassification$sp.b.distr$args)
    uncorrelated <- data.frame(se.as.uncor, sp.as.uncor, se.bs.uncor, 
        sp.bs.uncor)
    normalMoments <- function(t) {
        as.list(c(mean = mean(t), sd = sd(t)))
    }
    normalMomentList.se <- alply(uncorrelated[, c("se.as.uncor", 
        "se.bs.uncor")], 2, normalMoments)
    normalMomentList.sp <- alply(uncorrelated[, c("sp.as.uncor", 
        "sp.bs.uncor")], 2, normalMoments)
    normalize <- function(x) {
        (x - mean(x))/sd(x)
    }
    uncorrelated.normal <- apply(uncorrelated, 2, normalize)
    cor.matrix.se <- matrix(c(1, misclassification$se.cor, misclassification$se.cor, 
        1), nrow = 2)
    cor.matrix.sp <- matrix(c(1, misclassification$sp.cor, misclassification$sp.cor, 
        1), nrow = 2)
    if (misclassification$se.cor < 1) {
        correlated.normal.se <- qnorm(rcopula.gauss(iter, cor.matrix.se))
    }
    if (misclassification$sp.cor < 1) {
        correlated.normal.sp <- qnorm(rcopula.gauss(iter, cor.matrix.sp))
    }
    deNormalize <- function(x, sampleMean, sampleSd) {
        (x * sampleSd + sampleMean)
    }
    if (misclassification$se.cor < 1) {
        correlated.se <- correlated.normal.se
        for (j in 1:ncol(correlated.se)) {
            correlated.se[, j] <- deNormalize(correlated.normal.se[, 
                j], normalMomentList.se[[j]][[1]], normalMomentList.se[[j]][[2]])
        }
        correlated.se <- data.frame(correlated.se)
    }
    else {
        correlated.se <- uncorrelated[, c("se.as.uncor", "se.as.uncor")]
    }
    names(correlated.se) <- c("se.a", "se.b")
    if (misclassification$sp.cor < 1) {
        correlated.sp <- correlated.normal.sp
        for (j in 1:ncol(correlated.sp)) {
            correlated.sp[, j] <- deNormalize(correlated.normal.sp[, 
                j], normalMomentList.sp[[j]][[1]], normalMomentList.sp[[j]][[2]])
        }
        correlated.sp <- data.frame(correlated.sp)
    }
    else {
        correlated.sp <- uncorrelated[, c("sp.as.uncor", "sp.as.uncor")]
    }
    names(correlated.sp) <- c("sp.a", "sp.b")
    result <- data.frame(correlated.se, correlated.sp)
    if (iter.equaled.one) {
        result <- result[1, ]
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
