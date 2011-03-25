require(roxygen)

# Set working directory
directory <- "C:/Users/jthetzel/Research/pba/"
setwd(directory)
directory <- "/home/jthetzel/Research/pba/"
setwd(directory)


package.skeleton('pba', code_files='devel/R/pba.R', force=T, namespace=T, path=".")
roxygenize('pba', copy.package=T)




cat('"C:/Program Files/R/R-2.12.2/bin/x64/R" CMD build "C:/Users/jthetzel/Research/test.roxygen/"')
string <- '"C:/Program Files/R/R-2.12.2/bin/x64/R" CMD build "C:/Users/jthetzel/Research/test.roxygen/"'
system(string)

















#' Define bias parameters
#' 
#' @param variable A character vector naming the variable upon which the bias 
#' is effecting.
#' @param misclassification A list describing misclassification bias. The list 
#' consists of the following parameters:
#' \item {se.a.distr}{A pba.distr object defining the sensitivity among cases.}
#' \item {sp.a.distr}{A pba.distr object defining the specificity among cases.}
#' \item {se.b.distr}{A pba.distr object defining the sensitivity among 
#' non-cases.}
#' \item {sp.b.distr}{A pba.distr object defining the sensitivity among 
#' non-cases.}
#' \item {se.cor}{Correlation between sensitivity among cases and non-cases. 
#' Correlation of 1 indicates non-differential selection bias. Correlation of 
#' 0 indicates independent differential selection bias. Correlation less than 
#' 1 but greater than 0 indicates partial differential selection bias.}
#' \item {se.cor}{Correlation between specificity among cases and non-cases. 
#' Correlation of 1 indicates non-differential selection bias. Correlation of 0 
#' indicates independent differential selection bias. Correlation less than 1 
#' but greater than 0 indicates partial differential selection bias.}
#' @param selection A list describing selection bias. The list consists of the 
#' following parameters:
#' \item {s.a1.distr}{A pba.distr object defining selection among exposed cases.}
#' \item {s.a0.distr}{A pba.distr object defining selection among non-exposed 
#' cases.}
#' \item {s.b1.distr}{A pba.distr object defining selection among exposed 
#' non-cases.}
#' \item {s.b0.distr}{A pba.distr object defining selection among non-exposed 
#' non-cases.}
#' @param confounding A list containing one or more lists describing unmeasured 
#' confounding bias.  The confounder bias lists contain the following parameters:
#' \item {p1.distr}{A pba.distr object defining the probability of the 
#' unmeasured confounder among the exposed.}
#' \item {p0.distr}{A pba.distr object defining the probability of the 
#' unmeasured confounder among the non-exposed.}
#' \item {rr.distr}{(optional) A pba.distr object defining the relative risk 
#' association between the confounder and the outcome. rr.distr or or.distr must 
#' be defined. If both are defined, rr.distr is used instead of or.distr.
#' \item {or.distr}{(optional) A pba.distr object defining the odds ratio 
#' association between the confounder and the outcome. rr.distr or or.distr 
#' must be defined. If both are defined, rr.distr is used instead of or.distr.}
#' 
#' @return The pba function returns an object of class \code{pba.variables}.
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# Function to create a define bias of a variable
pbaVariable <- function(variable, # Character name of variable 
		misclassification = list(se.a.distr=NULL, sp.a.distr=NULL, 
				se.b.distr=NULL, sp.b.distr=NULL,
				se.cor=NULL, sp.cor=NULL),
		selection         = list(s.a1.distr=NULL, s.a0.distr=NULL, 
				s.b1.distr=NULL, s.b0.distr=NULL,
				s.1.cor=NULL, s.0.cor=NULL,
				or=NULL),
		confounding       = list(confounder=list(p1.distr=NULL, p0.distr=NULL, 
						rr.distr=NULL, rd.distr=NULL,
						name=NULL))
)
{											
	result <- list()
	result[[variable[[1]]]] <- list(variable=variable, 
			misclassification=misclassification, 
			selection=selection,
			confounding=confounding)
	
	class(result) <- "pba.variables"
	return(result)	
}