\name{plotBias}
\alias{plotBias}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotBias(pba, density = T, scales = "free", types = NULL, print = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pba}{
%%     ~~Describe \code{pba} here~~
}
  \item{density}{
%%     ~~Describe \code{density} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{types}{
%%     ~~Describe \code{types} here~~
}
  \item{print}{
%%     ~~Describe \code{print} here~~
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
function (pba, density = T, scales = "free", types = NULL, print = T) 
{
    if (is.null(types)) {
        types <- laply(pba$bias.tables, function(x) {
            names(x)
        })
    }
    types <- unique(types)
    data <- melt(pba$bias.tables)
    columns <- which(names(data) \%in\% c("variable", "value", 
        "L1", "L2", "L3"))
    data <- data[, columns]
    plots <- list()
    if ("misclassification" \%in\% types) {
        plots$misclassification <- pbaPlotMisclassification(data = data, 
            density = density, scales = scales)
    }
    if ("selection" \%in\% types) {
        plots$selection <- pbaPlotSelection(data = data, density = density, 
            scales = scales)
    }
    if ("confounding" \%in\% types) {
        plots$confounding.proportions <- pbaPlotConfoundingProportions(data = data, 
            density = density, scales = scales)
    }
    if ("confounding" \%in\% types) {
        plots$confounding.risks <- pbaPlotConfoundingRisks(data = data, 
            density = density, scales = scales)
    }
    if (print) {
        n <- length(types)
        vp <- list()
        j <- 1
        for (i in types) {
            vp[[i]] <- viewport(x = (j - 0.5)/n, width = 1/n, 
                just = c("centre", "centre"))
            j <- j + 1
        }
        grid.newpage()
        pushViewport(viewport())
        popViewport()
        for (i in names(plots)) {
            if (i \%in\% c("misclassification", "selection")) {
                pushViewport(vp[[i]])
                grid.draw(ggplotGrob(plots[[i]]))
                popViewport()
            }
            else if (i \%in\% c("confounding.proportions", "confounding.risks")) {
                pushViewport(vp[["confounding"]])
                if (i == "confounding.proportions") {
                  pushViewport(viewport(y = 1/3, height = 2/3, 
                    just = c("centre", "bottom")))
                  grid.draw(ggplotGrob(plots[[i]]))
                  popViewport()
                }
                if (i == "confounding.risks") {
                  pushViewport(viewport(y = 1/3, height = 1/3, 
                    just = c("centre", "top")))
                  grid.draw(ggplotGrob(plots[[i]]))
                  popViewport()
                }
                popViewport()
            }
        }
    }
    invisible(plots)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
