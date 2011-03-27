\name{pbaBiasListsTables}
\alias{pbaBiasListsTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaBiasListsTables(bias.lists)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bias.lists}{
%%     ~~Describe \code{bias.lists} here~~
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
function (bias.lists) 
{
    bias.tables <- list()
    for (i in names(bias.lists)) {
        misclassification.included <- ifelse(length(bias.lists[[i]][[1]]$misclassification) > 
            0, T, F)
        selection.included <- ifelse(length(bias.lists[[i]][[1]]$selection) > 
            0, T, F)
        confounding.included <- ifelse(length(bias.lists[[i]][[1]]$confounding) > 
            0, T, F)
        if (misclassification.included) {
            misclassification <- ldply(bias.lists[[i]], function(x) {
                x$misclassification
            })
            bias.tables[[i]]$misclassification <- misclassification
        }
        if (selection.included) {
            selection <- ldply(bias.lists[[i]], function(x) {
                x$selection
            })
            bias.tables[[i]]$selection <- selection
        }
        if (confounding.included) {
            confounding <- ldply(bias.lists[[i]], function(x) {
                ldply(x$confounding, function(y) {
                  y
                })
            })
            confounding <- split(confounding, f = confounding$.id)
            bias.tables[[i]]$confounding <- confounding
        }
    }
    return(bias.tables)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
