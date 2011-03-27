\name{print.pba.variables}
\alias{print.pba.variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.pba.variables(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    if (length(x) == 1) {
        variables <- paste(": ", names(x), sep = "")
    }
    else if (length(x) > 1) {
        variables <- paste("s: ", paste(names(x), collapse = ", "), 
            sep = "")
    }
    else {
        cat("\n\nNo bias parameters defined]n")
        break()
    }
    cat(paste("\nA pba.variables object describing bias parameters for the following\n\t\t\t\t\t\t\tvariable", 
        variables, "\n", sep = ""))
    for (i in names(x)) {
        cat(paste("\nVariable ", i, ":\n", sep = ""))
        if (!is.null(x[[i]]$misclassification[[1]])) {
            cat("\tMisclassification with:\n")
            cat(paste("\t\tSensitivity among cases distribution of", 
                x[[i]]$misclassification$se.a.distr$distr, "\n"))
            cat(paste("\t\tSensitivity among non-cases distribution of", 
                x[[i]]$misclassification$se.b.distr$distr, "\n"))
            cat(paste("\t\tSpecificity among cases distribution of", 
                x[[i]]$misclassification$sp.a.distr$distr, "\n"))
            cat(paste("\t\tSpecificity among non-cases distribution of", 
                x[[i]]$misclassification$sp.b.distr$distr, "\n"))
        }
        if (!is.null(x[[i]]$selection[[1]])) {
            cat("\tSelection with\n")
            cat(paste("\t\tSelection among exposed cases distribution of", 
                x[[i]]$selection$s.a1.distr$distr, "\n"))
            cat(paste("\t\tSelection among non-exposed cases distribution", 
                x[[i]]$selection$s.a0.distr$distr, "\n"))
            cat(paste("\t\tSelection among exposed non-cases distribution of", 
                x[[i]]$selection$s.b1distr$distr, "\n"))
            cat(paste("\t\tSelection among non-exposed non-cases distribution of", 
                x[[i]]$selection$s.b0.distr$distr, "\n"))
        }
        if (!is.null(x[[i]]$confounding[[1]][[1]])) {
            for (j in names(x[[i]]$confounding)) {
                cat(paste("\tUnmeasured confounding by", x[[i]]$confounding[[j]]$name, 
                  "with:\n"))
                cat(paste("\t\tProportion among exposed distribution of", 
                  x[[i]]$confounding[[j]]$p1.distr$distr, "\n"))
                cat(paste("\t\tProportion among non-exposed distribution of", 
                  x[[i]]$confounding[[j]]$p0.distr$distr, "\n"))
                if (!is.null(x[[1]]$confounding[[j]]$rr.distr)) {
                  cat(cat(paste("\t\tRelative risk of confounding distribution of", 
                    x[[i]]$confounding[[j]]$rr.distr$distr, "\n")))
                }
                if (!is.null(x[[1]]$confounding[[j]]$rd.distr)) {
                  cat(cat(paste("\t\tRisk difference of confounding distribution of", 
                    x[[i]]$confounding[[j]]$rr.distr$distr, "\n")))
                }
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
