\name{plotEstimates}
\alias{plotEstimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotEstimates(pba = NULL, data = NULL, density = T, exp = F, adjust = 1, binwidth = NULL, scales = "free", variables = NULL, print = T)
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
  \item{exp}{
%%     ~~Describe \code{exp} here~~
}
  \item{adjust}{
%%     ~~Describe \code{adjust} here~~
}
  \item{binwidth}{
%%     ~~Describe \code{binwidth} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{variables}{
%%     ~~Describe \code{variables} here~~
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
function (pba = NULL, data = NULL, density = T, exp = F, adjust = 1, 
    binwidth = NULL, scales = "free", variables = NULL, print = T) 
{
    if (!is.null(variables)) {
        data <- pba$coefficients.hat.random[variables]
    }
    if (is.null(variables)) {
        data <- pba$coefficients.hat.random
    }
    if (exp) {
        data <- llply(data, function(x) {
            exp(x)
        })
    }
    if (!density) {
        data <- melt(data)
        xlim <- quantile(data$value, c(0.01, 0.99))
        p1 <- ggplot(data, aes(x = value))
        p2 <- p1 + xlim(xlim) + facet_grid(L1 ~ ., scales = scales)
        plot <- p2 + geom_histogram(binwidth = binwidth)
    }
    if (density) {
        data <- ldply(data, function(x) {
            q.low <- quantile(x, 0.01)
            q.high <- quantile(x, 0.99)
            result <- density(x[x > q.low & x < q.high], adjust = adjust)
            lower <- as.numeric(quantile(x, 0.025))
            upper <- as.numeric(quantile(x, 0.975))
            data.frame(x = result$x, y = result$y, lower = lower, 
                upper = upper)
        })
        p1 <- ggplot(data, aes(x = x, y = y))
        p2 <- p1 + geom_line()
        ribbon <- geom_ribbon(data = subset(data, x >= lower & 
            x <= upper), aes(ymax = y), ymin = 0, alpha = 0.5)
        plot <- p2 + ribbon
        plot <- plot + scale_y_continuous("density") + scale_x_continuous("estimate") + 
            facet_grid(.id ~ ., scales = scales)
    }
    if (exp) {
        plot <- plot + scale_x_log(name = "estimate")
    }
    if (print) {
        print(plot)
    }
    invisible(plot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
