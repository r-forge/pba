\name{pbaIterateConfoundingSingle}
\alias{pbaIterateConfoundingSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbaIterateConfoundingSingle(exposure, outcome, model, confounding, name, iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exposure}{
%%     ~~Describe \code{exposure} here~~
}
  \item{outcome}{
%%     ~~Describe \code{outcome} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{confounding}{
%%     ~~Describe \code{confounding} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
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
function (exposure, outcome, model, confounding, name, iter) 
{
    models.updated <- list()
    for (i in 1:iter) {
        model.updated <- model
        data <- model.updated$data
        formula <- model.updated$formula
        rows.a1 <- which(data[, outcome] == 1 & data[, exposure] == 
            1)
        rows.a0 <- which(data[, outcome] == 1 & data[, exposure] == 
            0)
        rows.b1 <- which(data[, outcome] == 0 & data[, exposure] == 
            1)
        rows.b0 <- which(data[, outcome] == 0 & data[, exposure] == 
            0)
        a1. <- length(rows.a1)
        a0. <- length(rows.a0)
        b1. <- length(rows.b1)
        b0. <- length(rows.b0)
        expected <- pbaBackCalculateConfounding(a1. = a1., a0. = 0, 
            b1. = b1., b0. = b0., p1 = confounding$p1[i], p0 = confounding$p0[i], 
            rr = confounding$rr[i], rd = confounding$rd[i])
        prob.a1. <- expected$a11/a1.
        prob.a0. <- expected$a01/a0.
        prob.b1. <- expected$b11/b1.
        prob.b0. <- expected$b01/b0.
        confounder.a1 <- rbinom(a1., 1, prob.a1.)
        confounder.a0 <- rbinom(a0., 1, prob.a0.)
        confounder.b1 <- rbinom(b1., 1, prob.b1.)
        confounder.b0 <- rbinom(b0., 1, prob.b0.)
        confounder <- vector(length = nrow(data))
        confounder[rows.a1] <- confounder.a1
        confounder[rows.a0] <- confounder.a0
        confounder[rows.b1] <- confounder.b1
        confounder[rows.b0] <- confounder.b0
        data[[name]] <- confounder
        model.updated$data <- data
        model.updated$formula <- update.formula(model.updated$formula, 
            paste(deparse(model.updated$formula), " + ", name, 
                sep = ""))
        models.updated[[i]] <- model.updated
    }
    return(models.updated)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
