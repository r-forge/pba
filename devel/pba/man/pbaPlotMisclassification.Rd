\name{pbaPlotMisclassification}
\alias{pbaPlotMisclassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaPlotMisclassification(pba = NULL, data = NULL, density = T, parameters = c("se.a", "se.b", "sp.a", "sp.b"), title = "Misclassification", scales = "free", L2 = "misclassification")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pba}{
%%     ~~Describe \code{pba} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{density}{
%%     ~~Describe \code{density} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{L2}{
%%     ~~Describe \code{L2} here~~
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
function (pba = NULL, data = NULL, density = T, parameters = c("se.a", 
    "se.b", "sp.a", "sp.b"), title = "Misclassification", scales = "free", 
    L2 = "misclassification") 
{
    if (is.null(data)) {
        data <- melt(pba$bias.tables)
        columns <- which(names(data) \%in\% c("variable", "value", 
            "L1", "L2", "L3"))
        data <- data[, columns]
    }
    data <- subset(data, L2 == L2 & variable \%in\% parameters)
    facet <- facet_grid(variable ~ L1, scales = scales)
    plot1 <- ggplot(data, aes(x = value)) + facet
    plot1.density <- plot1 + geom_density()
    plot1.histogram <- plot1 + geom_histogram()
    if (density) {
        plot.result <- plot1.density
    }
    else {
        plot.result <- plot1.histogram
    }
    return(plot.result + opts(title = title))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
