# pba: Probabilistic Bias Analysis
# 
# Jeremy Thoms Hetzel
# jthetzel@gmail.com
###############################################################################

# Convenience functions: expit()
expit <- function(x)
{
	exp(x) / (1 + exp(x))
}

# Convenience functions: logit()
logit <- function(x)
{
	log(x / (1-x))
}



# Main function to perform PBA
pba <- function(model,
		pba.variables,
		iter = 1000,
		alpha = 0.05)
{
	# Create a data frame of bias estimates using pbaContigency function.
	# bias is a list of tables
	bias <- pbaMisclassification(model, pba.variables, iter)
	
	# Check if PPVs or NPVs are < 0 or > 1
	# If present, resample to create valid values
	bias <- pbaMisclassificationCheck(bias, model, pba.variables)
	
	# Predict unbiased dataset by iterating through data frame of bias estimates.
	# model.summaries is a list of model summaries
	model.summaries <- pbaIterate(model, bias$tables, iter)
	
	# Collect all parameter estimates from predicted model summaries
	coefficients.hat <- list()
	for (i in rownames(model.summaries[[1]]$coefficients))
	{
		coefficients.hat[[i]] <- t(sapply(model.summaries, function(x)
						{
							x$coefficients[i,]	
						}))
	}
	
	# Ajust expected coefficient estimates for random error
	coefficients.hat.random <- list()
	for (i in names(coefficients.hat))
	{
		coefficients.hat.random[[i]] <- coefficients.hat[[i]][,'Estimate'] + 
				rnorm(nrow(coefficients.hat[[i]]), 0, 1) * 
				coefficients.hat[[i]][,'Std. Error']
	}
	
	
	# Coefficients from original model
	coefficients.star <- summary(model)$coefficients
	
	# Create summary list of observed, bias adjusted, and bias and random
	# error adjusted
	summary <- pbaSummary(model=model,
			coefficients.star=coefficients.star,
			coefficients.hat=coefficients.hat,
			coefficients.hat.random=coefficients.hat.random,
			alpha=alpha)	
	
	# Summarize estimates with random and bias error
	estimates <- c()
	result <- list(bias=bias,
			coefficients.star=coefficients.star,
			coefficients.hat=coefficients.hat,
			coefficients.hat.random=coefficients.hat.random,
			summary=summary)
	class(result) <- c('pba', 'list')
	return(result)
	
}


# Main function to perform PBA
pba2 <- function(model,
		pba.variables,
		iter = 1000,
		alpha = 0.05)
{
	# Create a data frame of bias estimates using pbaContigency function.
	# bias is a list of tables
	bias.tables <- pbaBiasTables(pba.variables, iter)
	
	# Calculate predictive values for misclassification biases
	bias.tables <- pbaCalculatePredictiveValues(bias.tables = bias.tables,
			model = model, iter = iter)
	
	# Correct unrealistic misclassification values
	bias.tables <- pbaCorrectMisclassification(bias.tables = bias.tables,
			pbaVariables = pbaVariables, iter = iter)
	
	# Predict unbiased dataset by iterating through data frame of bias estimates.
	# model.summaries is a list of model summaries
	model.summaries <- pbaIterate(model, bias$tables, iter)
	
	# Collect all parameter estimates from predicted model summaries
	coefficients.hat <- list()
	for (i in rownames(model.summaries[[1]]$coefficients))
	{
		coefficients.hat[[i]] <- t(sapply(model.summaries, function(x)
						{
							x$coefficients[i,]	
						}))
	}
	
	# Ajust expected coefficient estimates for random error
	coefficients.hat.random <- list()
	for (i in names(coefficients.hat))
	{
		coefficients.hat.random[[i]] <- coefficients.hat[[i]][,'Estimate'] + 
				rnorm(nrow(coefficients.hat[[i]]), 0, 1) * 
				coefficients.hat[[i]][,'Std. Error']
	}
	
	
	# Coefficients from original model
	coefficients.star <- summary(model)$coefficients
	
	# Create summary list of observed, bias adjusted, and bias and random
	# error adjusted
	summary <- pbaSummary(model=model,
			coefficients.star=coefficients.star,
			coefficients.hat=coefficients.hat,
			coefficients.hat.random=coefficients.hat.random,
			alpha=alpha)	
	
	# Summarize estimates with random and bias error
	estimates <- c()
	result <- list(bias=bias,
			coefficients.star=coefficients.star,
			coefficients.hat=coefficients.hat,
			coefficients.hat.random=coefficients.hat.random,
			summary=summary)
	class(result) <- c('pba', 'list')
	return(result)
	
}

# Simulate classification of model data by Bernoulli trials
pbaIterate <- function(model, bias.tables, iter)
{
	# Create list to store glm summaries
	summaries <- list()
	
	# Iterations
	for (i in 1:iter)
	{
		# Create copy of model data to use during iterations
		data <- model$model
		
		# Iterate through bias tables
		for (j in names(bias.tables))
		{
			# Assign current table
			current.table <- bias.tables[[j]]
			
			# Find rows for a1s, a0s, b1s, and b0s
			rows.a1 <- which(model$model[,1]==1 & model$model[,j]==1)
			rows.a0 <- which(model$model[,1]==1 & model$model[,j]==0)
			rows.b1 <- which(model$model[,1]==0 & model$model[,j]==1)
			rows.b0 <- which(model$model[,1]==0 & model$model[,j]==0)
			
			
			# Simulate occurence of correct classification for a1s, a0s, b1s, and b0s
			correct.a1 <- rbinom(length(rows.a1), 1, current.table$ppv.a[i])
			correct.a0 <- rbinom(length(rows.a0), 1, current.table$npv.a[i])
			correct.b1 <- rbinom(length(rows.b1), 1, current.table$ppv.b[i])
			correct.b0 <- rbinom(length(rows.b0), 1, current.table$npv.b[i])
			
			# Change exposure if classification not correct
			data[rows.a1,j][correct.a1==0] <- as.numeric(!data[rows.a1,j][correct.a1==0])
			data[rows.a0,j][correct.a0==0] <- as.numeric(!data[rows.a0,j][correct.a0==0])
			data[rows.b1,j][correct.b1==0] <- as.numeric(!data[rows.b1,j][correct.b1==0])
			data[rows.b0,j][correct.b0==0] <- as.numeric(!data[rows.b0,j][correct.b0==0])
		}
		
		# Update model
		model.new <- update(model, data=data)
		
		# Save model summary
		summaries[[i]] <- summary(model.new)
	}
	
	return(summaries)
}


# Simulate classification of model data by Bernoulli trials
pbaIterate2 <- function(model, bias.tables, iter)
{
	# Create list to store glm summaries
	summaries <- list()
	
	# Store formula, incase it needs to be updated for additional unmeasured
	# confounders
	formula <- model$formula
	
	# Iterations
	for (i in 1:iter)
	{
		# Create copy of model data to use during iterations
		data <- model$model
		
		# Iterate through bias tables
		for (j in names(bias.tables))
		{
			# Assign current table
			current.table <- bias.tables[[j]]
			
			# Misclassification
			data <- pbaIterateMisclassification(exposure = j, model = model, 
					misclassification = current.table$misclassification, data = data,
					iter = iter)
			
			# Confouding
			for (k in names(current.table$confounding))
			{
				data <- pbaIterateConfounding(exposure = j, data = data,
						confounding = current.table$confounding[[k]],
						name = k, iter = iter)
				formula <- paste(formula, " + ", current.table$confounding[[k]]$name)			
			}
			
		}
			
			
		
		# Update model
		#! need to update call with additional confounders
		model.new <- update(model, data = data, formula = formula)
		
		# Save model summary
		summaries[[i]] <- summary(model.new)
	}
	
	return(summaries)
	
}


#
pbaIterateConfounding <- function(exposure, data, confounding, name, iter)
{
	# Find rows for a1s, a0s, b1s, and b0s
	rows.a1 <- which(data[,1]==1 & data[,exposure]==1)
	rows.a0 <- which(data[,1]==1 & data[,exposure]==0)
	rows.b1 <- which(data[,1]==0 & data[,exposure]==1)
	rows.b0 <- which(data[,1]==0 & data[,exposure]==0)
	
	# Calculate a1., a0., b1., b0.
	a1. <- length(rows.a1)
	a0. <- length(rows.a0)
	b1. <- length(rows.b1)
	b0. <- length(rows.b0)
	
	# Calculate expected confounding
	expected <- pbaBackCalculateCounfounding(a1. = a1., a0. = 0., b1. = b1., 
			b0. = b0., p1 = confounding$p1[i], p0 = confounding$p0[i], 
			rr = confounding$rr[i], rd = confounding$rd[i])
	
	# Calculate probabilities of confounding
	prob.a1. <- expected$a11 / a1.
	prob.a0. <- expected$a01 / a0.
	prob.b1. <- expected$b11 / b1.
	prob.b0. <- expected$b01 / b0.
	
	# Simulate occurence of correct classification for a1s, a0s, b1s, and b0s
	confounder.a1 <- rbinom(a1., 1, prob.a1.)
	confounder.a0 <- rbinom(a0., 1, prob.a0.)
	confounder.b1 <- rbinom(b1., 1, prob.b1.)
	confounder.b0 <- rbinom(b0., 1, prob.b0.)
	
	# Update data
	confounder <- vector(length = nrow(data))
	confounder[rows.a1] <- confounder.a1
	confounder[rows.a0] <- confounder.a0
	confounder[rows.b1] <- confounder.b1
	confounder[rows.b0] <- confounder.b0
	data[[name]] <- confounder
	
	
}



#
pbaIterateMisclassification <- function(exposure, model, misclassification, 
		data, iter)
{
	# Find rows for a1s, a0s, b1s, and b0s
	rows.a1 <- which(model$model[,1]==1 & model$model[,exposure]==1)
	rows.a0 <- which(model$model[,1]==1 & model$model[,exposure]==0)
	rows.b1 <- which(model$model[,1]==0 & model$model[,exposure]==1)
	rows.b0 <- which(model$model[,1]==0 & model$model[,exposure]==0)
	
	
	# Simulate occurence of correct classification for a1s, a0s, b1s, and b0s
	correct.a1 <- rbinom(length(rows.a1), 1, misclassification$ppv.a[i])
	correct.a0 <- rbinom(length(rows.a0), 1, misclassification$npv.a[i])
	correct.b1 <- rbinom(length(rows.b1), 1, misclassification$ppv.b[i])
	correct.b0 <- rbinom(length(rows.b0), 1, misclassification$npv.b[i])
	
	# Change exposure if classification not correct
	data[rows.a1,exposure][correct.a1==0] <- 
			as.numeric(!data[rows.a1,exposure][correct.a1==0])
	data[rows.a0,exposure][correct.a0==0] <- 
			as.numeric(!data[rows.a0,exposure][correct.a0==0])
	data[rows.b1,exposure][correct.b1==0] <- 
			as.numeric(!data[rows.b1,exposure][correct.b1==0])
	data[rows.b0,exposure][correct.b0==0] <- 
			as.numeric(!data[rows.b0,exposure][correct.b0==0])
	
	return(data)
}


#
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
#
#	  RETURNS Data frame with columns:
#			a1.hat
#			a0.hat
#			b1.hat
#			b0.hat
#			se.a
#			sp.a
#			se.b
#			sp.b
#			ppv.a
#			npv.a
#			ppv.b
#			npv.b	
{
	# Calculate false positive probability and false negative probability
	b.fp <- 1 - sp.b
	b.fn <- 1 - se.b
	a.fp <- 1 - sp.a
	
	# Calculate totoal cases and non-cases
	m1 <- a1.star + a0.star
	m0 <- b1.star + b0.star
	
	# Calculate expected counts
	b1.hat <- (sp.b * b1.star - b.fp*b0.star) / (se.b * sp.b - b.fn * b.fp)
	b0.hat <- m0 - b1.hat
	a1.hat <- (a1.star - a.fp * m1) / (se.a + sp.a - 1)
	a0.hat <- m1 - a1.hat
	
	# Calculate positive predictive value and negative predictive value
	ppv.a <- (se.a * a1.hat) / ((se.a * a1.hat) + ((1 - sp.a) * a0.hat))
	npv.a <- (sp.a * a0.hat) / ((sp.a * a0.hat) + ((1 - se.a) * a1.hat))
	ppv.b <- (se.b * b1.hat) / ((se.b * b1.hat) + ((1 - sp.b) * b0.hat))
	npv.b <- (sp.b * b0.hat) / ((sp.b * b0.hat) + ((1 - se.b) * b1.hat))
	
	# Return result
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



# Sample misclassification parameters
pbaSampleMisclassification <- function(misclassification, iter)
{
	# Create results list to store results
	results <- list()
	
	# Add number of iterations to 'n' argument in pba.variables object
	misclassification <- pbaAddIter(x=misclassification, iter=iter)
	
	# Sample the sensitivites and specificities
	se.as.uncor <- do.call(misclassification$se.a.distr$distr, 
			misclassification$se.a.distr$args)
	sp.as.uncor <- do.call(misclassification$sp.a.distr$distr, 
			misclassification$sp.a.distr$args)									
	se.bs.uncor <- do.call(misclassification$se.b.distr$distr, 
			misclassification$se.b.distr$args)
	sp.bs.uncor <- do.call(misclassification$sp.b.distr$distr, 
			misclassification$sp.b.distr$args)
	
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
	cor.matrix.se <- matrix(c(1, misclassification$se.cor,
					misclassification$se.cor, 1),
			nrow=2)
	cor.matrix.sp <- matrix(c(1, misclassification$sp.cor,
					misclassification$sp.cor, 1),
			nrow=2)
	
	# Generate correlated random draws with Gaussian copula
	# Skip if correlation equals 1
	if (misclassification$se.cor < 1) 
	{
		correlated.normal.se <- qnorm(rcopula.gauss(iter, cor.matrix.se))
	}
	if (misclassification$sp.cor < 1)
	{
		correlated.normal.sp <- qnorm(rcopula.gauss(iter, cor.matrix.sp))
	}
	
	# De-normalize the correlated random draws with the normal moments 
	# of the uncorrelated sensitivities and specificities
	deNormalize <- function(x, sampleMean, sampleSd) {
		(x * sampleSd + sampleMean)
	}	
	# First, sensitivities
	if (misclassification$se.cor < 1) 
	{
		correlated.se <- correlated.normal.se
		for (j in 1:ncol(correlated.se)) 
		{
			correlated.se[,j] <- deNormalize(correlated.normal.se[,j],
					normalMomentList.se[[j]][[1]],
					normalMomentList.se[[j]][[2]])
		}
		correlated.se <- data.frame(correlated.se)
	} else
	{
		correlated.se <- uncorrelated[,c('se.as.uncor', 'se.as.uncor')]
	}
	names(correlated.se) <- c('se.a', 'se.b')
	# Then, specificities
	if (misclassification$sp.cor < 1)
	{
		correlated.sp <- correlated.normal.sp
		for (j in 1:ncol(correlated.sp)) 
		{
			correlated.sp[,j] <- deNormalize(correlated.normal.sp[,j],
					normalMomentList.sp[[j]][[1]],
					normalMomentList.sp[[j]][[2]])
		}
		correlated.sp <- data.frame(correlated.sp)
	} else
	{
		correlated.sp <- uncorrelated[,c('sp.as.uncor', 'sp.as.uncor')]
	}
	names(correlated.sp) <- c('sp.a', 'sp.b')
	
	result <- data.frame(correlated.se, correlated.sp)
	
	# Return results list
	return(result)	
}

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
			# Add iter to n argument
			i <- pbaAddIter(i, iter)
			
			# Sample the proportion distributions
			p1 <- do.call(i$p1.distr$distr, i$p1.distr$args)
			p0 <- do.call(i$p1.distr$distr, i$p1.distr$args)
			
			# Correct if p1 is < or > 0 or 1, respectively
			p1 <- pbaCorrectConfounding(p = p1, distr = i$p1.distr$distr,
					args = i$p1.distr$args)
			p0 <- pbaCorrectConfounding(p = p0, distr = i$p0.distr$distr,
					args = i$p0.distr$args)
			
			# Sample rr if specified, otherwise sample rd
			if (!is.null(i$rr.distr))
			{
				rr <- do.call(i$rr.distr$distr, i$rr.distr$args)
				rd <- rep(NA, iter)
			} else
			{
				rd <- do.call(i$rd.distr$distr, i$rd.distr$args)
				rr <- rep(NA, iter)
			}
		}
		
		results[[i$name]] <- data.frame(p1, p0, rr, rd)	
	}

	return(results)
}






# Function to backcalculate confounding
pbaBackCalculateConfounding <- function(a1.=NULL, a0.=NULL, b1.=NULL, b0.=NULL, 
		p1=NULL, p0=NULL, rr=NULL, rd=NULL)
{
	# Calculate total exposed / non exposed, and confounded / non confounded
	n11 <- (a1. + b1.) * p1
	n01 <- (a0. + b0.) * p0
	n10 <- (a1. + b1.) - n11
	n00 <- (a0. + b0.) - n10
	n1.  <- a1. + b1.
	n0.  <- a0. + b0.
	
	if (!is.null(rr))
	{
		# Calculate expected counts by rr
		a11 <- (rr * n11 * a1.) / (rr * n11 + n1. - n11)
		a01 <- (rr * n01 * a0.) / (rr * n01 + n0. - n01)
		a10 <- a1. - a11
		a00 <- a0. - a01
		b11 <- n11 - a11
		b01 <- n01 - a01
		b10 <- b1. - b11
		b00 <- b0. - b01
	} else
	{
		# Calculate expected counts by rd
		a11 <- ((rd * n1. * (n1. - n11)) + (a1. * n11)) / n1.
		a01 <- ((rd * n0. * (n0. - n01)) + (a0. * n01)) / n0.
		a10 <- a1. - a11
		a00 <- a0. - a01
		b11 <- n11 - a11
		b01 <- n01 - a01
		b10 <- b1. - b11
		b00 <- b0. - b01	
	}
	
	# Return result
	result <- data.frame(a11, a01, a10, a00, b11, b01, b10, b00)
	return(result)
}







# Function to create a define bias of a variable
pbaVariable <- function(variable, # Character name of variable 
		misclassification = list(se.a.distr=NULL, sp.a.distr=NULL, 
														 se.b.distr=NULL, sp.b.distr=NULL,
														 se.cor=NULL, sp.cor=NULL),
		selection         = list(se.a.distr=NULL, sp.a.distr=NULL, 
													 	 se.b.distr=NULL, sp.b.distr=NULL,
														 se.cor=NULL, sp.cor=NULL,
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
	
	return(result)	
}



# Function to calculate sensitivity and specificity
pbaMisclassification <- function(model, pba.variables, iter)
# Returns list of tables with expected counts, sensitivity, specificity,
# positive predictive value, negative predictive value, observed odds ratio,
# and expected odds ratio
{
	# Create a results list to store results
	results <- list()
	
	# Create a replaced list to record number of rows that needed replacement due
	# to negative counts
	replaced <- list()
	
	# Adjust sensitivites and specificites for correlations
	# Correlated object is a list of tables
	correlated.tables <- pbaBiasCor(pba.variables, iter)
	
	# Iterate through correlated tables
	for (i in names(correlated.tables))
	{
		# Create observed contigency table from model
		table.star <- table(model$model[,1], model$model[,i])
		a1.star <- table.star[4]
		a0.star <- table.star[2]
		b1.star <- table.star[3]
		b0.star <- table.star[1]	
		
		# Create vector of counts to use with pbaBackCalculate
		a1.stars <- rep(a1.star, iter)
		a0.stars <- rep(a0.star, iter)
		b1.stars <- rep(b1.star, iter)
		b0.stars <- rep(b0.star, iter)
		
		# Back calculate expected counts, as well as
		# positive predicted values and negative predicted values
		result <- pbaBackCalculate(a1.star=a1.stars, 
				a0.star=a0.stars, 
				b1.star=b1.stars,
				b0.star=b0.stars,
				se.a = correlated.tables[[i]]$se.a,
				sp.a = correlated.tables[[i]]$sp.a,
				se.b = correlated.tables[[i]]$se.b,
				sp.b = correlated.tables[[i]]$sp.b)
		
		## Replace rows with counts less than 0	
		# Identify rows with counts less than 0
		rows <- which(apply(result[,c('a1.hat', 'a0.hat', 'b1.hat', 'b0.hat')], 1, 
						function(x)
						{
							any(x < 0 | is.na(x))
						}))
		
		# While rows exist with counts less than 0, replace the rows
		
		
		# Calculate the expected OR
		result$or.hat <- with(result, (a1.hat / a0.hat) / (b1.hat / b0.hat))
		
		# For comparison, calculate the observed OR
		result$or.star <- (a1.star / a0.star) / (b1.star / b0.star)
		
		# Store result in results list
		results[[i]] <- result
	}
	bias <- list(tables=results, replaced=replaced)
	
	return(bias)
}



# Remove values that result in PPV or NPV >1 or <0 and replace with new values
pbaMisclassificationCheck <- function(bias, model, pba.variables)
{
	names(bias)
	for(i in names(pba.variables))
	{
		replace <- T
		while (replace)
		{
			invalid <- apply(bias$tables[[i]][,c('se.a', 'sp.a', 'se.b', 'sp.b', 
									'ppv.a', 'npv.a', 'ppv.b', 'npv.b')], 1, 
					function(x)
					{
						any(c(is.na(x), x < 0, x > 1))
					})
			rows.replace <- which(invalid)
			
			iter.replace <- length(rows.replace)
			
			if (iter.replace > 0)
			{
				# Add plus one as a workaround for when
				# iter.replace = 1.  Otherwise, NAs result
				# because the pbaBiasCor  needs to calculate
				# standard deviation, and standard deviation of 1 is NA.
				bias.replace <- pbaMisclassification(model=model, 
						pba.variables=pba.variables[i], 
						iter=iter.replace+1)
				
				bias$tables[[i]][rows.replace,] <- bias.replace$tables[[i]][-(iter.replace + 1),]
				
				bias$replaced[[i]] <- bias$replaced[[i]] + bias.replace$replaced[[i]]
				
			} else
			{
				replace <- F
			}
		}
	}
	
	return(bias)
}



# Function to generate correlated random draws of 
# sensitivities and specificities.
# Adapted from JD Long: 
# http://www.cerebralmastication.com/2010/08/stochastic-simulation-with-copulas-in-r/
pbaBiasCor <- function(pba.variables, iter)
{
	# Create results list to store results
	results <- list()
	
	# Iterate through pba.variables 
	for (i in names(pba.variables))
	{
		# Assign this iteration's variable
		current.variable <- pba.variables[[i]]
		
		# Add number of iterations to 'n' argument in pba.variables object
		current.variable$misclassification$se.a.distr$args$n <- iter
		current.variable$misclassification$sp.a.distr$args$n <- iter
		current.variable$misclassification$se.b.distr$args$n <- iter
		current.variable$misclassification$sp.b.distr$args$n <- iter
		
		# Sample the sensitivites and specificities
		se.as.uncor <- do.call(current.variable$misclassification$se.a.distr$distr, 
				current.variable$misclassification$se.a.distr$args)
		sp.as.uncor <- do.call(current.variable$misclassification$sp.a.distr$distr, 
				current.variable$misclassification$sp.a.distr$args)									
		se.bs.uncor <- do.call(current.variable$misclassification$se.b.distr$distr, 
				current.variable$misclassification$se.b.distr$args)
		sp.bs.uncor <- do.call(current.variable$misclassification$sp.b.distr$distr, 
				current.variable$misclassification$sp.b.distr$args)
		
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
		cor.matrix.se <- matrix(c(1, current.variable$misclassification$se.cor,
						current.variable$misclassification$se.cor, 1),
				nrow=2)
		cor.matrix.sp <- matrix(c(1, current.variable$misclassification$sp.cor,
						current.variable$misclassification$sp.cor, 1),
				nrow=2)
		
		# Generate correlated random draws with Gaussian copula
		# Skip if correlation equals 1
		if (current.variable$misclassification$se.cor < 1) 
		{
			correlated.normal.se <- qnorm(rcopula.gauss(iter, cor.matrix.se))
		}
		if (current.variable$misclassification$sp.cor < 1)
		{
			correlated.normal.sp <- qnorm(rcopula.gauss(iter, cor.matrix.sp))
		}
		
		# De-normalize the correlated random draws with the normal moments 
		# of the uncorrelated sensitivities and specificities
		deNormalize <- function(x, sampleMean, sampleSd) {
			(x * sampleSd + sampleMean)
		}	
		# First, sensitivities
		if (current.variable$misclassification$se.cor < 1) 
		{
			correlated.se <- correlated.normal.se
			for (j in 1:ncol(correlated.se)) 
			{
				correlated.se[,j] <- deNormalize(correlated.normal.se[,j],
						normalMomentList.se[[j]][[1]],
						normalMomentList.se[[j]][[2]])
			}
			correlated.se <- data.frame(correlated.se)
		} else
		{
			correlated.se <- uncorrelated[,c('se.as.uncor', 'se.as.uncor')]
		}
		names(correlated.se) <- c('se.as', 'se.bs')
		# Then, specificities
		if (current.variable$misclassification$sp.cor < 1)
		{
			correlated.sp <- correlated.normal.sp
			for (j in 1:ncol(correlated.sp)) 
			{
				correlated.sp[,j] <- deNormalize(correlated.normal.sp[,j],
						normalMomentList.sp[[j]][[1]],
						normalMomentList.sp[[j]][[2]])
			}
			correlated.sp <- data.frame(correlated.sp)
		} else
		{
			correlated.sp <- uncorrelated[,c('sp.as.uncor', 'sp.as.uncor')]
		}
		names(correlated.sp) <- c('sp.a', 'sp.b')
		
		result <- data.frame(uncorrelated, correlated.se, correlated.sp)
		
		# Assign result to results list
		results[[i]] <- result
	}	
	
	# Return results list
	names(results) <- names(pba.variables)
	return(results)
}


# Method to print pbaVariables object
#print.pbaVariables <- function(pba.variables)
#{
#	print(paste('pbaVariables object for variable', variable, sep=''))
#	print(paste('Sensitivity of '))
#	
#}

# Method to print pba object
print.pba <- function(pba)
{
	return(pba$summary)
}


# Create summary list of observed, bias adjusted, and bias and random
# error adjusted
pbaSummary <- function(model, coefficients.star, coefficients.hat, 
		coefficients.hat.random, alpha)
{
	# Create summary list to store results
	summary <- list()
	
	# Summarize estimate and confidience limits from observed model
	summary$star <- data.frame(estimate=coefficients.star[,'Estimate'],
			confint(model, level=1 - alpha),
			p=coefficients.star[,'Pr(>|z|)'])
	names(summary$star) <- c('estimate', 
			paste(100 * alpha / 2,'%', sep=''),
			paste(100 * (1 - alpha / 2),'%', sep=''),
			'Pr(>|z|)')
	
	# Summarize expected estimate and simulation limits after bias adjustment
	summary$hat <- t(sapply(coefficients.hat, function(x)
					{
						quantile(x[,'Estimate'], c(0.5, alpha / 2, 1 - alpha / 2))	
					}))
	summary$hat <- data.frame(summary$hat)
	names(summary$hat) <- c('estimate', 
			paste(100 * alpha / 2,'%', sep=''),
			paste(100 * (1 - alpha / 2),'%', sep=''))
	# Simulated p value from simulations with estimate < 0	
	summary$hat$p <- as.vector(sapply(coefficients.hat, function(x)
					{
						length(x[,'Estimate'][x[,'Estimate'] < 0]) / length(x[,'Estimate'])
					}))
	summary$hat$p <- ifelse(summary$hat$p > 0.5, 
			1 - summary$hat$p, summary$hat$p)
	
	# Summarize expected estimate and simulation limits after bias and
	# random error adjustment
	summary$hat.random <- t(sapply(coefficients.hat.random, function(x)
					{
						quantile(x, c(0.5, alpha / 2, 1 - alpha / 2))	
					}))
	summary$hat.random <- data.frame(summary$hat.random)
	names(summary$hat.random) <- c('estimate', 
			paste(100 * alpha / 2,'%', sep=''),
			paste(100 * (1 - alpha / 2),'%', sep=''))
	# Simulated p value from simulations with estimate < 0
	summary$hat.random$p <- as.vector(sapply(coefficients.hat.random, function(x)
					{
						length(x[x<0]) / length(x)
					}))
	summary$hat.random$p <- ifelse(summary$hat.random$p > 0.5, 
			1 - summary$hat.random$p, summary$hat.random$p)
	
	# Return summary
	return(summary)
}




# Remove values that result in PPV or NPV >1 or <0 and replace with new values
pbaBiasCorCheck <- function(result, pba.variables)
{
	for(i in names(pba.variables))
	{
		replace <- T
		while (replace)
		{
			invalid <- apply(result[[i]][,c('ppv.a', 'npv.a', 'ppv.b', 'npv.b')], 1, 
					function(x)
					{
						any(c(is.na(x), x < 0, x > 1))
					})
			rows.replace <- which(invalid)
			
			iter.replace <- length(rows.replace)
			
			if (iter.replace > 0)
			{
				result.replace <- pbaBiasCor(pba.variables=pba.variables[i], 
						iter=iter.replace+1) # Add plus one as a workaround for when
				# iter.replace = 1.  Otherwise, NAs result
				# because the pbaBiasCor  needs to calculate
				# standard deviation.
				result[[i]][rows.replace,] <- result.replace[-(iter.replace+1),]
			} else
			{
				replace <- F
			}
		}
	}
	
	return(result)
}


# Plot distribution of simulated estimates
pbaPlotEstimates <- function(pba, density=T, exp=F, adjust=1, binwidth=NULL,
		scales='free', variables=NULL)
{
	if (!is.null(variables))
	{
		data <- pba$coefficients.hat.random[variables]
	}
	
	if (is.null(variables))
	{
		data <- pba$coefficients.hat.random
	}
	
	# Apply exp() transformation
	if (exp)
	{
		data <- llply(data, function(x)
				{
					exp(x)
				})
	}
	
	# Histogram
	if (!density)
	{
		data <- melt(data)
		
		xlim <- quantile(data$value, c(0.01, 0.99)) # Trim outliers
		p1 <- ggplot(data, aes(x=value))
		p2 <- p1 + xlim(xlim) +	facet_grid(L1~., scales=scales)
		plot <- p2 + geom_histogram(binwidth=binwidth)
	}
	
	# Density
	if (density)
	{
		data <- ldply(data, function(x)
				{
					q.low <- quantile(x, 0.01) # Trim off lower outliers
					q.high <- quantile(x, 0.99) # Trim off upper outliers
					result <- density(x[x > q.low & x < q.high], adjust=adjust) # Density
					lower <- as.numeric(quantile(x, 0.025))
					upper <- as.numeric(quantile(x, 0.975))
					data.frame(x=result$x, y=result$y, lower=lower, 
							upper=upper)
				})
		
		p1 <- ggplot(data, aes(x=x, y=y))
		p2 <- p1 + geom_line()
		ribbon <- geom_ribbon(data=subset(data, x >= lower & x <= upper), 
				aes(ymax=y), ymin=0, alpha=0.5)
		plot <- p2 + ribbon	
		
		# Scales and facets
		plot <- plot + scale_y_continuous("density") +
				scale_x_continuous("estimate") +
				facet_grid(.id~., scales=scales)
	}
	
	# Transform x axis if exp
	if (exp)
	{
		plot <- plot + scale_x_log(name="estimate") 
	}					 					 
	
	# Return plot
	return(plot)
}

# Plot sensitivities and specificities
# Biases
pbaPlotBias <- function(pba, density=T)
{
	# Include only relevant columns
	data <- lapply(pba$bias$tables, function(x)
			{
				x[c('se.a', 'se.b', 'sp.a', 'sp.b')]
			})
	
	# Melt data into long-form
	data <- melt(data)
	
	# Density
	p1 <- ggplot(data, aes(x=value))
	p2 <- p1 + geom_density() + facet_grid(variable~L1)
	
	# Histogram
	p3 <- p1 + geom_histogram() + facet_grid(variable~L1)
	
	if (density)
	{
		return(p2)
	} else
	{
		return(p3)
	}
}



pbaAddIter <- function(x, iter)
{
	for (i in names(x))
	{
		if (is.list(x[[i]]))
		{
			x[[i]]$args$n <- iter
		}
	}
	
	return(x)
}




# Funciton to generate bias tables from pbaVariables
pbaBiasTables <- function(pbaVariables, iter)
{
	results <- list()
	
	for (i in pbaVariables)
	{
		misclassification <- pbaSampleMisclassification(
			misclassification = i$misclassification, iter = iter)
		selection <- list()
		confounding <- pbaSampleConfounding(
			confounding = i$confounding, iter = iter)
		results[[i$variable]] <- list(misclassification = misclassification,
																selection = selection,
																confounding = confounding)
	}
	
	return(results)
}


# Function to calculate PPV and NPV for misclassification
pbaCalculatePredictiveValues <- function(bias.tables, model, iter)
{
	for (i in names(bias.tables))
	{
		result <- pbaCalculatePredictiveValuesInternal(exposure = i,
				misclassification = bias.tables[[i]]$misclassification, model = model,
				iter = iter)
		# Bind result to bias table
		bias.tables[[i]]$misclassification <- cbind(bias.tables[[i]]$misclassification,
				result[,c('a1.hat','a0.hat','b1.hat','b0.hat','ppv.a','npv.a','ppv.b',
				'npv.b')])
	}	

	return(bias.tables)
}


#
pbaCalculatePredictiveValuesInternal <- function(exposure, misclassification, 
		model, iter)
{
	# Create observed contigency table from model
	table.star <- table(model$model[,1], model$model[,exposure])
	a1.star <- table.star[4]
	a0.star <- table.star[2]
	b1.star <- table.star[3]
	b0.star <- table.star[1]
	
	# Create vector of counts to use with pbaBackCalculate
	a1.stars <- rep(a1.star, iter)
	a0.stars <- rep(a0.star, iter)
	b1.stars <- rep(b1.star, iter)
	b0.stars <- rep(b0.star, iter)
	
	# Back calculate expected counts, as well as
	# positive predicted values and negative predicted values
	result <- pbaBackCalculate(a1.star=a1.stars, 
			a0.star=a0.stars, 
			b1.star=b1.stars,
			b0.star=b0.stars,
			se.a = misclassification$se.a,
			sp.a = misclassification$sp.a,
			se.b = misclassification$se.b,
			sp.b = misclassification$sp.b)
	
	return(result)
}





# Remove and resample if values are less than 0 or greater than 1
pbaCorrectConfounding <- function(p, distr, args)
{
	rows <- which(p < 0 | p > 1)
	while (length(rows) > 0)
	{
		args$n <- length(rows)
		p[rows] <- do.call(distr, args)
		rows <- which(p < 0 | p > 1)
	}
	
	return(p)
}

# Remove and resample if misclassification
pbaCorrectMisclassification <- function(bias.tables, pbaVariables, model, iter)
{
	for (i in names(pbaVariables))
	{
		# Find rows needing replacement
		rows <- which(bias.tables[[i]]$misclassification$a1.hat < 0 |	
						bias.tables[[i]]$misclassification$a0.hat < 0 |
						bias.tables[[i]]$misclassification$b1.hat < 0 |
						bias.tables[[i]]$misclassification$b0.hat < 0 |
						bias.tables[[i]]$misclassification$ppv.a < 0 |
						bias.tables[[i]]$misclassification$ppv.a > 1 |
						bias.tables[[i]]$misclassification$npv.a < 0 |
						bias.tables[[i]]$misclassification$npv.a > 1 |
						bias.tables[[i]]$misclassification$ppv.b < 0 |
						bias.tables[[i]]$misclassification$ppv.b > 1 |
						bias.tables[[i]]$misclassification$npv.b < 0 |
						bias.tables[[i]]$misclassification$npv.b > 1 |
						bias.tables[[i]]$misclassification$se.a < 0 |
						bias.tables[[i]]$misclassification$se.a > 1 |
						bias.tables[[i]]$misclassification$se.b < 0 |
						bias.tables[[i]]$misclassification$se.b > 1 |
						bias.tables[[i]]$misclassification$sp.a < 0 |
						bias.tables[[i]]$misclassification$sp.a > 1 |
						bias.tables[[i]]$misclassification$sp.b < 0 |
						bias.tables[[i]]$misclassification$sp.b > 1 )
		
		while (length(rows) > 0)
		{
			# Resample sensitivities and specificities for rows needing replacement
			bias.tables[[i]]$misclassification[rows,c('se.a','se.b','sp.a','sp.b')] <- 
					pbaSampleMisclassification(
					misclassification = pbaVariables[[i]]$misclassification,
					iter = length(rows))
	
			# Calculate predictive values of replaced rows
			result <- pbaCalculatePredictiveValuesInternal(exposure = i,
					misclassification = bias.tables[[i]]$misclassification[rows,], 
					model = model, iter = length(rows))
			
			# Replace predictive values
			bias.tables[[i]]$misclassification[rows,c('a1.hat','a0.hat','b1.hat',
							'b0.hat','ppv.a','npv.a','ppv.b','npv.b')] <- 
					result[,c('a1.hat','a0.hat','b1.hat','b0.hat','ppv.a','npv.a','ppv.b',
									'npv.b')]
			
			# Find rows needing replacement
			rows <- which(bias.tables[[i]]$misclassification$a1.hat < 0 |	
							bias.tables[[i]]$misclassification$a0.hat < 0 |
							bias.tables[[i]]$misclassification$b1.hat < 0 |
							bias.tables[[i]]$misclassification$b0.hat < 0 |
							bias.tables[[i]]$misclassification$ppv.a < 0 |
							bias.tables[[i]]$misclassification$ppv.a > 1 |
							bias.tables[[i]]$misclassification$npv.a < 0 |
							bias.tables[[i]]$misclassification$npv.a > 1 |
							bias.tables[[i]]$misclassification$ppv.b < 0 |
							bias.tables[[i]]$misclassification$ppv.b > 1 |
							bias.tables[[i]]$misclassification$npv.b < 0 |
							bias.tables[[i]]$misclassification$npv.b > 1 |
							bias.tables[[i]]$misclassification$se.a < 0 |
							bias.tables[[i]]$misclassification$se.a > 1 |
							bias.tables[[i]]$misclassification$se.b < 0 |
							bias.tables[[i]]$misclassification$se.b > 1 |
							bias.tables[[i]]$misclassification$sp.a < 0 |
							bias.tables[[i]]$misclassification$sp.a > 1 |
							bias.tables[[i]]$misclassification$sp.b < 0 |
							bias.tables[[i]]$misclassification$sp.b > 1 )
		}
	}
	
	return(bias.tables)
}
