# Sample misclassification parameters
pbaSampleMisclassification <- function(misclassification, iter)
{	
	# Create results list to store results
	results <- list()
	
	# Specify sigma, which is the lower triangle of a correlation matrix
	if (is.matrix(misclassification$sigma))
	{
		sigma <- lower.tri(misclassification$sigma)
	} else if (is.vector(misclassification$sigma))
	{
		sigma <- misclassification$sigma	
	} else
	{
		sigma <- c(misclassification$se.cor, 0, 0, 0, 0, misclassification$sp.cor)
	}
	
	# Generate correlated samples
	correlated <- pbaCorrelate(n = iter, pbaDistrs = misclassification, 
			sigma = sigma)
	
	
	# Return results data frame
	return(data.frame(correlated))	
}





pbaCorrelate <- function (n, pbaDistrs, sigma)
{
	# Extract distributions
	distr.list <- lapply(pbaDistrs[which(lapply(pbaDistrs, class)=="pbaDistr")], 
			function(x) 
				{
					x$distr
				})
	# Extraxt arguments
	args.list <- lapply(pbaDistrs[which(lapply(pbaDistrs, class)=="pbaDistr")], 
			function(x) 
			{
				x$args
			})
	
	# Create normal copula to construct correlations
	copula <- normalCopula(param = sigma, dim = length(distr.list), dispstr = "un")
	
	# Create multivariate distribution
	correlated <- mvdc(copula = copula, margins = paste(distr.list), paramMargins = args.list)
	
	out <- rmvdc(correlated, n = n)
	colnames(out) <- gsub(".distr", "", 
			names(which(lapply(pbaDistrs, class)=="pbaDistr")))
	
	return(out)
}


	
	
	