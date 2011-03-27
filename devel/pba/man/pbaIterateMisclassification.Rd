\name{pbaIterateMisclassification}
\alias{pbaIterateMisclassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaIterateMisclassification(model, bias.tables, iter)
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
    for (i in 1:iter) {
        model.updated <- model
        data <- model.updated$data
        for (j in names(bias.tables)) {
            current.table <- bias.tables[[j]]
            rows.a1 <- which(model$model[, 1] == 1 & model$model[, 
                j] == 1)
            rows.a0 <- which(model$model[, 1] == 1 & model$model[, 
                j] == 0)
            rows.b1 <- which(model$model[, 1] == 0 & model$model[, 
                j] == 1)
            rows.b0 <- which(model$model[, 1] == 0 & model$model[, 
                j] == 0)
            correct.a1 <- rbinom(length(rows.a1), 1, current.table$misclassification$ppv.a[i])
            correct.a0 <- rbinom(length(rows.a0), 1, current.table$misclassification$npv.a[i])
            correct.b1 <- rbinom(length(rows.b1), 1, current.table$misclassification$ppv.b[i])
            correct.b0 <- rbinom(length(rows.b0), 1, current.table$misclassification$npv.b[i])
            data[rows.a1, j][correct.a1 == 0] <- as.numeric(!data[rows.a1, 
                j][correct.a1 == 0])
            data[rows.a0, j][correct.a0 == 0] <- as.numeric(!data[rows.a0, 
                j][correct.a0 == 0])
            data[rows.b1, j][correct.b1 == 0] <- as.numeric(!data[rows.b1, 
                j][correct.b1 == 0])
            data[rows.b0, j][correct.b0 == 0] <- as.numeric(!data[rows.b0, 
                j][correct.b0 == 0])
        }
        model.updated$data <- data
        models.updated[[i]] <- model.updated
    }
    return(models.updated)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
