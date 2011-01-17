# TODO: Add comment
# 
# Author: jthetzel
###############################################################################


# Method from Rothman and colleagues, differing normal, updated 2
exp$se.cor <- .75
exp$sp.cor <- 1
pbaBias <- exp
iter <- 1000
	# Add number of iterations to 'n' argument in pbaBias object
	pbaBias$se.a.distr$args$n <- iter
	pbaBias$sp.a.distr$args$n <- iter
	pbaBias$se.b.distr$args$n <- iter
	pbaBias$sp.b.distr$args$n <- iter
	
	# Sample the sensitivites and specificities
	se.as.uncor <- do.call(pbaBias$se.a.distr$distr, pbaBias$se.a.distr$args)
	sp.as.uncor <- do.call(pbaBias$sp.a.distr$distr, pbaBias$sp.a.distr$args)									
	se.bs.uncor <- do.call(pbaBias$se.b.distr$distr, pbaBias$se.b.distr$args)
	sp.bs.uncor <- do.call(pbaBias$sp.b.distr$distr, pbaBias$sp.b.distr$args)
	
	# Combine sampled uncorrelated sensitivities and specificities 
  # into a data frame
	uncorrelated <- data.frame(se.as.uncor, sp.as.uncor, se.bs.uncor, sp.bs.uncor)
	
	# Record the normal moments in a list
	normalMoments <- function(t) {
		as.list(c(mean=mean(t), sd=sd(t)))
	}
	normalMomentList.se <- alply(uncorrelated[,c('se.as.uncor','se.bs.uncor')], 
			2, normalMoments)
	normalMomentList.sp <- alply(uncorrelated[,c('sp.as.uncor','sp.bs.uncor')], 
			2, normalMoments)
	
	# Normalize the uncorrelated sensitivities and specificities
	normalize <- function(x) {
		(x - mean(x)) / sd(x)
	}
	uncorrelated.normal <- apply(uncorrelated, 2, normalize)
	
	# Specify correlation matrix
	cor.matrix.se <- matrix(c(1, pbaBias$se.cor,
									 					pbaBias$se.cor, 1),
														nrow=2)
	cor.matrix.sp <- matrix(c(1, pbaBias$sp.cor,
														pbaBias$sp.cor, 1),
														nrow=2)
	
	# Generate correlated random draws with Gaussian copula
	# Skip if correlation equals 1
	if (pbaBias$se.cor < 1) 
	{
		correlated.normal.se <- qnorm(rcopula.gauss(iter, cor.matrix.se))
	}
	if (pbaBias$sp.cor < 1)
	{
		correlated.normal.sp <- qnorm(rcopula.gauss(iter, cor.matrix.sp))
	}

	# De-normalize the correlated random draws with the normal moments 
  # of the uncorrelated sensitivities and specificities
	deNormalize <- function(x, sampleMean, sampleSd) {
		(x * sampleSd + sampleMean)
	}	
	# First, sensitivities
	if (pbaBias$se.cor < 1) 
	{
		correlated.se <- correlated.normal.se
		for (i in 1:ncol(correlated.se)) 
		{
			correlated.se[,i] <- deNormalize(correlated.normal.se[,i],
					normalMomentList.se[[i]][[1]],
					normalMomentList.se[[i]][[2]])
		}
		correlated.se <- data.frame(correlated.se)
	} else
	{
		correlated.se <- uncorrelated[,c('se.as.uncor', 'se.as.uncor')]
	}
		names(correlated.se) <- c('se.as', 'se.bs')
		# Second, specificities
		if (pbaBias$sp.cor < 1)
		{
		correlated.sp <- correlated.normal.sp
		for (i in 1:ncol(correlated.sp)) 
		{
			correlated.sp[,i] <- deNormalize(correlated.normal.sp[,i],
					normalMomentList.sp[[i]][[1]],
					normalMomentList.sp[[i]][[2]])
		}
		correlated.sp <- data.frame(correlated.sp)
	} else
	{
		correlated.sp <- uncorrelated[,c('sp.as.uncor', 'sp.as.uncor')]
	}
	names(correlated.sp) <- c('sp.as', 'sp.bs')
	
	result <- data.frame(uncorrelated, correlated.se, correlated.sp)
	
	# Remove rows where sensitivity or specificity are <0 or >1
	result <- subset(result, !se.as < 0 & !se.as > 1 &
					!se.bs < 0 & !se.bs > 1 &
					!sp.as < 0 & !sp.as > 1 &
					!sp.bs < 0 & !sp.bs > 1)
	
	# Warn if rows removed
	if(nrow(result) != iter) { warning(paste(iter - nrow(result),
						' iterations removed due to sensitivity or specificity values
								less than 0 or greater than 1.'))}
	
	return(result)
}

exp$se.cor <- 1
exp$sp.cor <- .8
test <- pbaBiasCor(exp, 1000)
cor(test[,-c(1:4)])

# try ggplot
require(ggplot2)

meltDraws <-melt(test[,-c(1:4)])
meltDraws$dataSource <- "simulated"
meltData <- melt(test[,c(1:4)])
meltData$dataSource <- "original"
plotData <- rbind(meltData, meltDraws)

qplot(value, colour=dataSource, data=plotData, geom="density")+ facet_wrap(~variable)

