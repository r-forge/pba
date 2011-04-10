# Sample confounding paramaters
pbaSampleConfounding <- function(confounding, iter)
{
	# Create a results list to store results
	results <- list()
	
	# Iterate through variables
	for (i in confounding)
	{
		if (!is.null(i$p1.distr))
		{
			# Specify sigma, which is the lower triangle of a correlation matrix
			if (is.matrix(i$sigma))
			{
				sigma <- lower.tri(i$sigma)
			} else if (is.vector(i$sigma))
			{
				sigma <- i$sigma	
			} else if (!is.null(i$p.cor))
			{
				sigma <- c(i$p.cor)
			} else
			{
				sigma <- 0
			}
			
			# Generate correlated samples
			correlated <- pbaCorrelate(n = iter, pbaDistrs = i[c('p1.distr', 'p0.distr')], 
					sigma = sigma)
						
			# Correct if p1 is < or > 0 or 1, respectively
			correlated <- pbaCorrectProbability(p = correlated, pbaDistr = i,
					sigma = sigma)
			
			# Sample rr if specified, otherwise sample rd
			if (!is.null(i$rr.distr))
			{
				args <- i$rr.distr$args
				args$n <- iter
				rr <- do.call(paste("r", i$rr.distr$distr, sep=""), args = args)
				rd <- rep(NA, iter)
			} else
			{
				args <- i$rd.distr$args
				args$n <- iter
				rd <- do.call(paste("r", i$rd.distr$distr, sep=""), args = args)
				rr <- rep(NA, iter)
			}
		}
		
		results[[i$name]] <- data.frame(correlated, rr, rd)	
	}
	
	return(results)
}
