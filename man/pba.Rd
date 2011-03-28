\name{pba}
\alias{pba}
\title{Main function to perform a probabilistic bias analysis...}
\usage{pba(model, pba.variables, iter=1000)}
\description{Main function to perform a probabilistic bias analysis}
\value{The pba function returns an object of class \code{pba}.
An object of class \code{pba} is a list containing at least the following 
components:
\item{call}{The matched call.}
\item{bias.tables}{A list of values of the bias parameters sampled at each 
iteration.}
\item{coefficients.hat}{A list of the coefficients at each iteration after 
adjusting for bias.}
\item{coefficients.hat.random}{A list of the coefficients at each iteration 
after adjusting for bias and random error.}
\item{coefficients.star}{Named vector of coefficients from the original glm 
or lm object, before adjustment for bias.}
\item{model}{The model frame from the original glm or lm object.}
\item{iter}{Number of iterations.}
\item{pba.variables}{The pba.variables object describing the bias 
parameters.}
\item{model.summaries}{List of model summaries at each iteration.}
\item{time.elapsed}{Time elapsed during evaluation of the pba function.}}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\arguments{\item{model}{A glm or lm object.}
\item{pba.variables}{One or more pba.variables objects which define bias 
parameters.}
\item{iter}{Number of iterations to perform.}}
