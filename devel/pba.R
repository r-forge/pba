# Probabilistic Bias Analysis
# 
# Jeremy T. Hetzel
# jthetzel@gmail.com
###############################################################################

# General functions
expit <- function(x)
{
	exp(x) / (1 + exp(x))
}

logit <- function(x)
{
	log(x / (1-x))
}



# Main function to perform PBA.
pba <- function(model,
								pbaBias,
								iter = 1000,
								alpha = 0.05
								)
{
	# Create a data frame of bias estimates using pbaContigency function
	bias.table <- pbaMisclassification(model, pbaBias, iter)
	
	# Predict unbiased dataset by iterating through data frame of bias estimates
	model.summaries <- pbaIterate(model, bias.table$ppv.a, bias.table$npv.a, 
																bias.table$ppv.b, bias.table$npv.b, iter)
																									
	# Collect all parameter estimates from predicted model summaries
	coefficients.hat <- list()
	for (i in rownames(model.summaries[[1]]$coefficients))
	{
		coefficients.hat[[i]] <- t(sapply(model.summaries, function(x)
		{
		 	x$coefficients[i,]	
		}))
	}
	
	# Coefficients from original model
	coefficients.star <- summary(model)$coefficients
	
	summary <- t(sapply(coefficients.hat, function(x)
					{
						quantile(x[,'Estimate'], c(0.5, alpha / 2, 1 - alpha / 2))	
					}))
	
	# Summarize estimates with random and bias error
	estimates <- c()
	result <- list(bias.table=bias.table,
								 coefficients.hat=coefficients.hat,
								 coefficients.star=coefficients.star,
								 summary=summary)

}


# Simulate classification of model data by Bernoulli trials
pbaIterate <- function(model, ppv.a, npv.a, ppv.b, npv.b, iter)
{
	# Create list to store glm summaries
	summaries <- list()
	
	# Iterations
	for (i in 1:iter)
	{
		# Create copy of model data to use during iterations
		data <- model$model
		
		# Find rows for a1s, a0s, b1s, and b0s
		rows.a1 <- which(model$model[,1]==1 & model$model[,2]==1)
		rows.a0 <- which(model$model[,1]==1 & model$model[,2]==0)
		rows.b1 <- which(model$model[,1]==0 & model$model[,2]==1)
		rows.b0 <- which(model$model[,1]==0 & model$model[,2]==0)
		
		# Simulate occurence of correct classification for a1s, a0s, b1s, and b0s
		correct.a1 <- rbinom(length(rows.a1), 1, ppv.a[i])
		correct.a0 <- rbinom(length(rows.a0), 1, npv.a[i])
		correct.b1 <- rbinom(length(rows.b1), 1, ppv.b[i])
		correct.b0 <- rbinom(length(rows.b0), 1, npv.b[i])
		
		# Change exposure if classification not correct
		# print(data[rows.b0, 2][correct.b0 == 0])
		data[rows.a1,2][correct.a1==0] <- as.numeric(!data[rows.a1,2][correct.a1==0])
		data[rows.a0,2][correct.a0==0] <- as.numeric(!data[rows.a0,2][correct.a0==0])
		data[rows.b1,2][correct.b1==0] <- as.numeric(!data[rows.b1,2][correct.b1==0])
		data[rows.b0,2][correct.b0==0] <- as.numeric(!data[rows.b0,2][correct.b0==0])
		
		# Update model
		model.new <- update(model, data=data)
		
		# Save model summary
		summaries[[i]] <- summary(model.new)
	}
	
	return(summaries)
}

# Function to perform PBA on contigency table
pbaContigency <- function(
		a1.star,
		a0.star,
		b1.star,
		b0.star,
		se.a.distr,
		sp.a.distr,
		se.b.distr,
		sp.b.distr,
		iter)
{
	se.a.distr$args$n <- iter
	sp.a.distr$args$n <- iter
	se.b.distr$args$n <- iter
	sp.b.distr$args$n <- iter
	
	se.as <- do.call(se.a.distr$distr, se.a.distr$args)
	sp.as <- do.call(sp.a.distr$distr, sp.a.distr$args)									
	se.bs <- do.call(se.b.distr$distr, se.b.distr$args)
	sp.bs <- do.call(sp.b.distr$distr, sp.b.distr$args)
	
	a1.stars <- rep(a1.star, iter)
	a0.stars <- rep(a0.star, iter)
	b1.stars <- rep(b1.star, iter)
	b0.stars <- rep(b0.star, iter)
	

	result <- pbaBackCalculate(a1.star=a1.stars, 
															 a0.star=a0.stars, 
														   b1.star=b1.stars,
														   b0.star=b0.stars,
								             	 se.a = se.as,
															 sp.a = sp.as,
															 se.b = se.bs,
															 sp.b = sp.bs)
												
	result$or.star <- (a1.star / a0.star) / (b1.star / b0.star)
	result$or.hat <- with(result, (a1.hat / a0.hat) / (b1.hat / b0.hat))
	
	return(result)
}




pbaBackCalculate <- function(a1.star, a0.star, b1.star, b0.star, 
														 se.a, sp.a, se.b, sp.b)
# Function to back calculate the true number of exposed and 
# unexposed subjects given the observed exposed, 
# observed unexposed, sensitivity, and specificity.
# Adapted from equations [19-7] and [19-8] in Rothman KJ et al.'s 
# Modern Epidemiology (2008).
#   ARGUMENTS:
#   	a1.star = observed exposed cases
#   	a0.star = observed non-exposed cases
#  		b1.star = observed exposed non-cases (or controls)
#  	  b0.star = observed non-exposed non-cases (or controls)
#			se.a    = sensitivity among cases
#			sp.a    = specificity among cases
#			se.b    = sensitivity among non-cases (or controls)
#			sp.b    = specificity among non-cases (or controls)
{
	b.fp <- 1 - sp.b
	b.fn <- 1 - se.b
	a.fp <- 1 - sp.a
	
	m1 <- a1.star + a0.star
	m0 <- b1.star + b0.star
	
	b1.hat <- (sp.b * b1.star - b.fp*b0.star) / (se.b * sp.b - b.fn * b.fp)
	b0.hat <- m0 - b1.hat
	a1.hat <- (a1.star - a.fp * m1) / (se.a + sp.a - 1)
	a0.hat <- m1 - a1.hat
	
	ppv.a <- (se.a * a1.hat) / ((se.a * a1.hat) + ((1 - sp.a) * a0.hat))
	npv.a <- (sp.a * a0.hat) / ((sp.a * a0.hat) + ((1 - se.a) * a1.hat))
	ppv.b <- (se.b * b1.hat) / ((se.b * b1.hat) + ((1 - sp.b) * b0.hat))
	npv.b <- (sp.b * b0.hat) / ((sp.b * b0.hat) + ((1 - se.b) * b1.hat))
	
	######################################
	############!!!!!!!!##################
	######################################
	# Set PPVs and NPVs greater than 1 to 1
	ppv.a[ppv.a > 1] <- 1
	npv.a[npv.a > 1] <- 1
	ppv.b[ppv.b > 1] <- 1
	npv.b[npv.b > 1] <- 1

	######################################
	############!!!!!!!!##################
	######################################
	# Set PPVs and NPVs less than 0 to 0
	ppv.a[ppv.a < 0] <- 0
	npv.a[npv.a < 0] <- 0
	ppv.b[ppv.b < 0] <- 0
	npv.b[npv.b < 0] <- 0
	
	result <- data.frame(a1.hat, a0.hat, b1.hat, b0.hat, se.a, sp.a, se.b, sp.b, 
											 ppv.a, npv.a, ppv.b, npv.b)
	return(result)
}



	
# Function to create a pba distribution object
pbaDistr <- function(distr, args)
# Character name of random generation
# distribution function. E.g. 'rnorm'.
# List of arguments passed to ditribution function.
# E.g. args=list(mean=0.85, sd=0.05) for the 
# rnorm function.  'n' should not be specified, as
# it will be automatically determined in subsequent
# PBA functions.
{
	result <- list(distr=distr, args=args)
}

# Function to create a define bias of a variable
pbaBias <- function(variable, # Character name of variable 
										se.a.distr, sp.a.distr, 
										se.b.distr, sp.b.distr,
										se.cor, sp.cor)
{
	result <- list(variable=variable, 
								 se.a.distr=se.a.distr, sp.a.distr=sp.a.distr, 
								 se.b.distr=se.b.distr, sp.b.distr=sp.b.distr, 
								 se.cor=se.cor, sp.cor=sp.cor)
	return(result)	
}



# Function to calculate sensitivity and specificity
pbaMisclassification <- function(model,
																 pbaBias,
																 iter)
{
	# Adjust sensitivites and specificites for correlations
	correlated <- pbaBiasCor(pbaBias, iter)
	
	# Create observed contigency table from model
	table.star <- table(model$model[,1], model$model[,pbaBias$variable])
	a1.star <- table.star[4]
	a0.star <- table.star[2]
	b1.star <- table.star[3]
	b0.star <- table.star[1]	
	
	# Create vector of counts to use with pbaBackCalculate
	a1.stars <- rep(a1.star, iter)
	a0.stars <- rep(a0.star, iter)
	b1.stars <- rep(b1.star, iter)
	b0.stars <- rep(b0.star, iter)
	
	# Back calculate expected counts
	result <- pbaBackCalculate(a1.star=a1.stars, 
			a0.star=a0.stars, 
			b1.star=b1.stars,
			b0.star=b0.stars,
			se.a = correlated$se.as,
			sp.a = correlated$sp.as,
			se.b = correlated$se.bs,
			sp.b = correlated$sp.bs)
	
	# Calculate the expected OR
	result$or.hat <- with(result, (a1.hat / a0.hat) / (b1.hat / b0.hat))
	
	# For comparison, calculate the observed OR
	result$or.star <- (a1.star / a0.star) / (b1.star / b0.star)
	
	return(result)
}

# Function to generate correlated random draws of 
# sensitivities and specificities.
# Adapted from JD Long: 
# http://www.cerebralmastication.com/2010/08/stochastic-simulation-with-copulas-in-r/
require('QRMlib')
require('plyr')
pbaBiasCor <- function(pbaBias, iter)
{
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
	
	######################################
	############!!!!!!!!##################
	######################################
	# Set sensitivities and specificities less than 0 to 0
	result[result<0] <- 0
	
	######################################
	############!!!!!!!!##################
	######################################
	# Set sensitivities and specificities greater than 1 to 1
	result[result>1] <- 1

	# Remove rows where sensitivity or specificity are <0 or >1
	#result <- subset(result, !se.as < 0 & !se.as > 1 &
	#				!se.bs < 0 & !se.bs > 1 &
	#				!sp.as < 0 & !sp.as > 1 &
	#				!sp.bs < 0 & !sp.bs > 1)
	
	# Warn if rows removed
	if(nrow(result) != iter) { warning(paste(iter - nrow(result),
						' iterations removed due to sensitivity or specificity values
								less than 0 or greater than 1.'))}
	
	return(result)
}



# Method to print pbaBias object
print.pbaBias <- function(pbaBias)
{
	print(paste('pbaBias object for variable', variable, sep=''))
	print(paste('Sensitivity of '))
	
}