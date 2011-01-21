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
#	  RETURNS Data frame with columsns:
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

# Function to create a define bias of a variable
pbaVariable <- function(variable, # Character name of variable 
		se.a.distr=NULL, sp.a.distr=NULL, 
		se.b.distr=NULL, sp.b.distr=NULL,
		se.cor=NULL, sp.cor=NULL)
{
	misclassification <- list(se.a.distr=se.a.distr, sp.a.distr=sp.a.distr, 
														se.b.distr=se.b.distr, sp.b.distr=sp.b.distr, 
														se.cor=se.cor, sp.cor=sp.cor)
												
	result <- list()
	result[[variable[[1]]]] <- list(variable=variable, misclassification=misclassification)
			
	return(result)	
}



# Function to calculate sensitivity and specificity
pbaMisclassification <- function(model,
																 pba.variables,
																 iter)
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
															 se.a = correlated.tables[[i]]$se.as,
															 sp.a = correlated.tables[[i]]$sp.as,
															 se.b = correlated.tables[[i]]$se.bs,
															 sp.b = correlated.tables[[i]]$sp.bs)
													 
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
		names(correlated.sp) <- c('sp.as', 'sp.bs')
		
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
			rows.replace <- which(invalide)
			
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
