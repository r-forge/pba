\name{pbaIterateSelection}
\alias{pbaIterateSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaIterateSelection(model, model.original = NULL, bias.tables, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{model.original}{
%%     ~~Describe \code{model.original} here~~
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
function (model, model.original = NULL, bias.tables, iter) 
{
    models.updated <- list()
    for (i in 1:iter) {
        model.updated <- model
        data <- model.updated$data
        rows <- nrow(data)
        selection <- llply(bias.tables, function(x) {
            x$selection
        })
        probs <- list()
        for (j in names(selection)) {
            if (length(selection[[j]]) > 0) {
                prob <- vector(length = rows)
                rows.a1 <- which(data[, 1] == 1 & data[, j] == 
                  1)
                rows.a0 <- which(data[, 1] == 1 & data[, j] == 
                  0)
                rows.b1 <- which(data[, 1] == 0 & data[, j] == 
                  1)
                rows.b0 <- which(data[, 1] == 0 & data[, j] == 
                  0)
                prob[rows.a1] <- selection[[j]]$s.a1
                prob[rows.a0] <- selection[[j]]$s.a0
                prob[rows.b1] <- selection[[j]]$s.b1
                prob[rows.b0] <- selection[[j]]$s.b0
                probs[[j]] <- prob
            }
        }
        not.selected.prob <- 1 - Reduce("*", probs)
        if (length(not.selected.prob) > 0) {
            not.selected <- rbinom(rows, 1, not.selected.prob)
            model.updated$data <- rbind(data, data[which(not.selected == 
                1), ])
            models.updated[[i]] <- model.updated
        }
        else {
            models.updated[[i]] <- model
        }
    }
    return(models.updated)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
