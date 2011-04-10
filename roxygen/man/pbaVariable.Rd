\name{pbaVariable}
\alias{pbaVariable}
\title{Define bias parameters...}
\usage{pbaVariable(variable, misclassification, selection, confounding)pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
pbaVariable(variable, misclassification, selection, confounding)
}
\description{Define bias parameters}
\value{The pba function returns an object of class \code{pba.variables}.\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


\code{pbaVariable}: The pba function returns an object of class \code{pba.variables}.


}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\arguments{\item{variable}{A character vector naming the variable upon which the bias 
is effecting.}
\item{misclassification}{A list describing misclassification bias. The list 
consists of the following parameters:
\item{se.a.distr}{A pba.distr object defining the sensitivity among cases.}
\item{sp.a.distr}{A pba.distr object defining the specificity among cases.}
\item{se.b.distr}{A pba.distr object defining the sensitivity among 
non-cases.}
\item{sp.b.distr}{A pba.distr object defining the sensitivity among 
non-cases.}
\item{se.cor}{Correlation between sensitivity among cases and non-cases. 
Correlation of 1 indicates non-differential selection bias. Correlation of 
0 indicates independent differential selection bias. Correlation less than 
1 but greater than 0 indicates partial differential selection bias.}
\item{se.cor}{Correlation between specificity among cases and non-cases. 
Correlation of 1 indicates non-differential selection bias. Correlation of 0 
indicates independent differential selection bias. Correlation less than 1 
but greater than 0 indicates partial differential selection bias.}}
\item{selection}{A list describing selection bias. The list consists of the 
following parameters:
\item{s.a1.distr}{A pba.distr object defining selection among exposed cases.}
\item{s.a0.distr}{A pba.distr object defining selection among non-exposed 
cases.}
\item{s.b1.distr}{A pba.distr object defining selection among exposed 
non-cases.}
\item{s.b0.distr}{A pba.distr object defining selection among non-exposed 
non-cases.}}
\item{confounding}{A list containing one or more lists describing unmeasured 
confounding bias.  The confounder bias lists contain the following parameters:
\item{p1.distr}{A pba.distr object defining the probability of the 
unmeasured confounder among the exposed.}
\item{p0.distr}{A pba.distr object defining the probability of the 
unmeasured confounder among the non-exposed.}
\item{rr.distr}{(optional) A pba.distr object defining the relative risk 
association between the confounder and the outcome. rr.distr or or.distr must 
be defined. If both are defined, rr.distr is used instead of or.distr.}
\item{or.distr}{(optional) A pba.distr object defining the odds ratio 
association between the confounder and the outcome. rr.distr or or.distr 
must be defined. If both are defined, rr.distr is used instead of or.distr.}}}
\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

\alias{pbaVariable}

