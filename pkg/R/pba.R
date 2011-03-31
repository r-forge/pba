###############################################################################
# pba: Probabilistic Bias Analysis                                            #
# version 0.2-1                                                               #
# Jeremy Thoms Hetzel                                                         #
# jthetzel@gmail.com                                                          #
###############################################################################


#' Performs bias analysis with probabilistic methods
#' 
#' \tabular{ll}{
#' Package: \tab pba\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2-1\cr
#' Date: \tab 2011-03-22\cr
#' Title: \tab Probabilistic Bias Analysis\cr
#' Author@@R: \tab person("Jeremy", "Hetzel", "Thoms",email = "jthetzel@@gmail.com")\cr
#' Author: \tab Jeremy T. Hetzel <jthetzel@@gmail.com>\cr
#' Maintainer: \tab Jeremy T. Hetzel \email{jthetzel@@gmail.com}\cr
#' Depends: \tab R (>= 2.12.0), plyr, ggplot2\cr
#' Description: \tab The pba package performs probabilistic bias analyses, 
#' with methods are adapted from Lash TL, Fox MP, and Fink AK's Applying 
#' Quantitative Bias Analysis to Epidemiological Data (2010). Currently,
#' only methods for adjustment of misclassification bias are included.\cr
#' License: \tab GPL-3\cr
#' URL: \tab http://pba.r-forge.r-project.org/\cr
#' BugReports: \tab https://r-forge.r-project.org/tracker/?func=add&group_id=988&atid=3913\cr		
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
NULL

#' pba: Probabilistic Bias Analysis
#' 
#' The pba package performs probabilistic bias analyses, with methods are 
#' adapted from Lash TL, Fox MP, and Fink AK's Applying Quantitative Bias 
#' Analysis to Epidemiological Data (2010). Currently, only methods for 
#' adjustment of misclassification bias are included.
#' 
#' @name pba-package
#' @docType package
#' @import plyr
#' @import reshape
#' @import ggplot2
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}
#' 
NULL



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


#' Main function to perform a probabilistic bias analysis
#' 
#' @param model A glm or lm object.
#' @param pba.variables One or more pba.variables objects which define bias 
#' parameters.
#' @param iter Number of iterations to perform.
#' 
#' @return The pba function returns an object of class \code{pba}.
#' An object of class \code{pba} is a list containing at least the following 
#' components:
#' \item{call}{The matched call.}
#' \item{bias.tables}{A list of values of the bias parameters sampled at each 
#' iteration.}
#' \item{coefficients.hat}{A list of the coefficients at each iteration after 
#' adjusting for bias.}
#' \item{coefficients.hat.random}{A list of the coefficients at each iteration 
#' after adjusting for bias and random error.}
#' \item{coefficients.star}{Named vector of coefficients from the original glm 
#' or lm object, before adjustment for bias.}
#' \item{model}{The model frame from the original glm or lm object.}
#' \item{iter}{Number of iterations.}
#' \item{pba.variables}{The pba.variables object describing the bias 
#' parameters.}
#' \item{model.summaries}{List of model summaries at each iteration.}
#' \item{time.elapsed}{Time elapsed during evaluation of the pba function.}
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# Main function to perform PBA
pba <- function(model,
		pba.variables,
		iter = 1000,
		progress = 100)
{
	# Record start time
	time.start <- Sys.time()
	
	call <- match.call()
	
	# Create summary of original model
	model.summary <- summary(model)
	
	# Create list to store updated models and bias tables
	model.summaries <- list()
	bias.lists <- list()
	
	for (i in 1:iter)
	{
		# Print iteration to follow progress
		if (is.null(progess) | is.na(progress))
		{
			if (i %% progress == 0) print(i)
		}
		
		# Sample bias parameters
		bias.tables <- pbaBiasTables(pba.variables, iter = 1)
		
		# Calculate predictive values for misclassification biases
		bias.tables <- pbaCalculatePredictiveValues(bias.tables = bias.tables,
				model = model, iter = 1)
		
		# Create copy of model to update
		model.updated <- model
		
		## Misclassification
		# Correct unrealistic misclassification values
		bias.tables <- pbaCorrectMisclassification(bias.tables = bias.tables,
				pba.variables = pba.variables, model = model, iter = 1)		
		
		# Adjust data for misclassification bias
		model.updated <- pbaIterateMisclassification(model = model.updated, 
				bias.tables = bias.tables, iter = 1)[[1]]
		
		## Selection bias
		# Adjust for selection bias
		model.updated <- pbaIterateSelection(model = model.updated,
				bias.tables = bias.tables, iter =1)[[1]]
		
		## Unmeasured confounding
		# Adjust data for misclassification
		model.updated <- pbaIterateConfounding(model = model.updated, 
				bias.tables = bias.tables, iter = 1)[[1]]
		
		## Summarize models
		# Update model
		model.updated <- update(model.updated, data = model.updated$data, 
				formula = model.updated$formula)
		
		
		# Save result and bias table
		model.summaries[[i]] <- summary(model.updated)['coefficients']
		for (j in names(pba.variables))
		{
			bias.lists[[j]][[i]] <- bias.tables[[j]]
		}
		
	}
	
	
	# Convert bias lists into bias tables
	bias.tables <- pbaBiasListsTables(bias.lists = bias.lists)
	
	# Collect all parameter estimates from model summaries
	coefficients.hat <- list()
	for (i in rownames(model.summaries[[1]]$coefficients))
	{
		coefficients.hat[[i]] <- t(sapply(model.summaries, function(x)
						{
							x$coefficients[i,]	
						}))
	}
	
	# Adjust expected coefficient estimates for random error
	coefficients.hat.random <- list()
	for (i in names(coefficients.hat))
	{
		coefficients.hat.random[[i]] <- coefficients.hat[[i]][,'Estimate'] + 
				rnorm(nrow(coefficients.hat[[i]]), 0, 1) * 
				rep(model.summary$coefficients[,'Std. Error'], nrow(coefficients.hat[[i]]))
	}
	
	# Coefficients from original model
	coefficients.star <- summary(model)$coefficients
	
	# Record end time
	time.end <- Sys.time()
	time.elapsed <- time.end - time.start
	
	results <- list(
			call = call,
			bias.tables = bias.tables,
			coefficients.hat = coefficients.hat,
			coefficients.hat.random = coefficients.hat.random,
			coefficients.star = coefficients.star,
			bias.lists = bias.lists,
			model = model,
			iter = iter,
			pba.variables = pba.variables,
			model.summaries = model.summaries,
			time.elapsed = time.elapsed)
	
	class(results) <- "pba"
	return(results)
}



# Simulate classification of model data by Bernoulli trials
pbaIterateMisclassification <- function(model, bias.tables, iter)
{
	# Create list to store updated models
	models.updated <- list()
	
	# Iterations
	for (i in 1:iter)
	{
		# Create copy of model to use during iterations
		model.updated <- model
		data <- model.updated$data
		
		# Iterate variables
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
			correct.a1 <- rbinom(length(rows.a1), 1, 
					current.table$misclassification$ppv.a[i])
			correct.a0 <- rbinom(length(rows.a0), 1, 
					current.table$misclassification$npv.a[i])
			correct.b1 <- rbinom(length(rows.b1), 1, 
					current.table$misclassification$ppv.b[i])
			correct.b0 <- rbinom(length(rows.b0), 1, 
					current.table$misclassification$npv.b[i])

			# Change exposure if classification not correct				
			data[rows.a1,j][correct.a1==0] <- as.numeric(!data[rows.a1,j][correct.a1==0])
			data[rows.a0,j][correct.a0==0] <- as.numeric(!data[rows.a0,j][correct.a0==0])
			data[rows.b1,j][correct.b1==0] <- as.numeric(!data[rows.b1,j][correct.b1==0])
			data[rows.b0,j][correct.b0==0] <- as.numeric(!data[rows.b0,j][correct.b0==0])
		}
		
		# Save updated model
		model.updated$data <- data
		models.updated[[i]] <- model.updated
	}
	
	return(models.updated)
}


# Remove and resample if misclassification
pbaCorrectMisclassification <- function(bias.tables, pba.variables, model, iter)
{
	for (i in names(pba.variables))
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
							misclassification = pba.variables[[i]]$misclassification,
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







#
pbaIterateConfoundingSingle <- function(exposure, outcome, model, confounding, 
		name, iter)
{
	# Create list to updated models
	models.updated <- list()
	
	for (i in 1:iter)
	{
		# Create copy of model data to use during iterations
		model.updated <- model
		data <- model.updated$data
		formula <- model.updated$formula
		
		# Find rows for a1s, a0s, b1s, and b0s
		rows.a1 <- which(data[,outcome]==1 & data[,exposure]==1)
		rows.a0 <- which(data[,outcome]==1 & data[,exposure]==0)
		rows.b1 <- which(data[,outcome]==0 & data[,exposure]==1)
		rows.b0 <- which(data[,outcome]==0 & data[,exposure]==0)
		
		# Calculate a1., a0., b1., b0.
		a1. <- length(rows.a1)
		a0. <- length(rows.a0)
		b1. <- length(rows.b1)
		b0. <- length(rows.b0)
		
		# Calculate expected confounding
		expected <- pbaBackCalculateConfounding(a1. = a1., a0. = 0., b1. = b1., 
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
		
		# Update model
		model.updated$data <- data
		model.updated$formula <- update.formula(model.updated$formula, 
				paste(deparse(model.updated$formula), ' + ', name, sep = ''))
		
		# Save updated data
		models.updated[[i]] <- model.updated
	}
	
	return(models.updated)
}


# 
pbaIterateConfounding <- function(model, bias.tables, iter)
{
	# Create list to store updated models
	models.updated <- list()
	
	# Make a copy of the model to update
	model.updated <- model
	
	# Specify name of outcome variable
	outcome <- names(model$model)[1]
	
	for (i in 1:iter)
	{	
		# Iterate through exposure variables
		for (i in names(bias.tables))
		{
			# Continue only if confounding specified in bias tables
			if (length(bias.tables[[i]]$confounding) > 0)
			{
				# Iterate through unmeasured confounders
				for (j in names(bias.tables[[i]]$confounding))
				{
					exposure <- i
					confounding <- bias.tables[[i]]$confounding[[j]]
					model.updated <- pbaIterateConfoundingSingle(exposure = exposure,
							outcome = outcome, model = model.updated, confounding = confounding,
							name = j, iter = 1)[[1]]			
				}
			}
		}
		
		# Store updated model
		models.updated[[i]] <- model.updated
	}
	
	return(models.updated)
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



# Convert bias lists into bias tables
pbaBiasListsTables <- function(bias.lists)
{
	bias.tables <- list()
	
	# Iterate through variables
	for (i in names(bias.lists))
	{
		# Check is misclassification, selection, confounding exist
		misclassification.included <- ifelse(
				length(bias.lists[[i]][[1]]$misclassification) > 0, T, F)
		selection.included <- ifelse(
				length(bias.lists[[i]][[1]]$selection) > 0, T, F)
		confounding.included <- ifelse(
				length(bias.lists[[i]][[1]]$confounding) > 0, T, F)
		
		# Misclassification
		if (misclassification.included)
		{
			misclassification <- ldply(bias.lists[[i]], function(x)
					{
						x$misclassification
					})
			bias.tables[[i]]$misclassification <- misclassification
		}
		
		# Selection
		if (selection.included)
		{
			selection <- ldply(bias.lists[[i]], function(x)
					{
						x$selection
					})
			bias.tables[[i]]$selection <- selection
		}
		
		# Confounding
		if (confounding.included)
		{
			confounding <- ldply(bias.lists[[i]], function(x)
					{
						ldply(x$confounding, function(y)
								{
									y
								})
					})
			
			confounding <- split(confounding, f = confounding$.id)
			bias.tables[[i]]$confounding <- confounding
		}
	}
	
	return(bias.tables)
}








# Funciton to generate bias tables from pba.variables object
pbaBiasTables <- function(pba.variables, iter)
{
	results <- list()
	
	for (i in pba.variables)
	{
		if (!is.null(i$misclassification$se.a.distr))
		{
			misclassification <- pbaSampleMisclassification(
					misclassification = i$misclassification, iter = iter)
		} else 
		{
			misclassification = list()
		}
		
		if (!is.null(i$selection$s.a1.distr))
		{
			selection <- pbaSampleSelection(selection=i$selection, iter=iter)
		} else
		{
			selection <- list()
		}
		
		if (!is.null(i$confounding[[1]]$p1.distr))
		{
			confounding <- pbaSampleConfounding(
					confounding = i$confounding, iter = iter)
		} else
		{
			confounding <- list()
		}		
		
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
		if (length(bias.tables[[i]]$misclassification) > 0)
		{
			result <- pbaCalculatePredictiveValuesInternal(exposure = i,
					misclassification = bias.tables[[i]]$misclassification, model = model,
					iter = iter)
			# Bind result to bias table
			bias.tables[[i]]$misclassification <- cbind(bias.tables[[i]]$misclassification,
					result[,c('a1.hat','a0.hat','b1.hat','b0.hat','ppv.a','npv.a','ppv.b',
									'npv.b')])
		}
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




# Sample selection paramaters
pbaSampleSelection <- function(selection, iter)
{
	# Create a results list to store results
	results <- c()
	
	if (!is.null(selection$s.a1.distr))
	{
		# Add iter to n argument
		selection <- pbaAddIter(selection, iter)
		
		# Sample the proportion distributions
		s.a1 <- do.call(selection$s.a1.distr$distr, selection$s.a1.distr$args)
		s.a0 <- do.call(selection$s.a0.distr$distr, selection$s.a0.distr$args)
		s.b1 <- do.call(selection$s.b1.distr$distr, selection$s.b1.distr$args)
		s.b0 <- do.call(selection$s.b0.distr$distr, selection$s.b0.distr$args)
		
		# Correct if p1 is < or > 0 or 1, respectively
		s.a1 <- pbaCorrectConfounding(p = s.a1, distr = selection$s.a1.distr$distr,
				args = selection$s.a1.distr$args)
		s.a0 <- pbaCorrectConfounding(p = s.a0, distr = selection$s.a0.distr$distr,
				args = selection$s.a0.distr$args)
		s.b1 <- pbaCorrectConfounding(p = s.b1, distr = selection$s.b1.distr$distr,
				args = selection$s.b1.distr$args)
		s.b0 <- pbaCorrectConfounding(p = s.b0, distr = selection$s.b0.distr$distr,
				args = selection$s.b0.distr$args)
		
	}
	result <- data.frame(s.a1, s.a0, s.b1, s.b0)	
	
	
	return(result)
}




# Sample misclassification parameters
pbaSampleMisclassification <- function(misclassification, iter)
{
	# Function will return NAs if iter == 1, I think do to the way
	# the correlations are imposed.  As a temporary fix, if iter ==
	# 1, then iter is changed to 2, and later on the first row of the 
	# result is returned.
	iter.equaled.one <- if(iter == 1) T else F
	iter <- if(iter == 1) 2 else iter
	
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
	
	# If iter originally equaled 1, then only return first row.
	# A second row was added to avoid the returning of NAs
	if (iter.equaled.one)
	{
		result <- result[1,]
	}
	
	# Return results data frame
	return(result)	
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


#' @param pba A \code{pba} pbject.
#' @param transformation (optional) A character string naming a function to 
#' apply to the coefficients. If using a function that transforms the 
#' coefficients to a multiplicative scale, be sure to specify the \code{scale} 
#' argument accordingly.
#' @param scale A character string of either \code{"additive"} or 
#' \code{"multiplicative"}. If \code{"additive"}, precision will be calculated 
#' the upper confidience level minus the lower confidence level. If 
#' \code{"multiplicative"}, the precision will be calculated as the upper 
#' confidence level divided by the lower confidence interval. Precision is 
#' calculated after any transformations.
#' @param alpha Alpha level used for determining confidence limits.
#' @param ... (optional) Additional arguments to pass to the 
#' \code{transformation} function, if specified.
#' 
#' @S3method print pba
#' @rdname pba
#' @export
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# Create summary list of observed, bias adjusted, and bias and random
# error adjusted
summary.pba <- function(pba, transformation=NULL, scale="additive", alpha=0.05, ...)
{
	
	model <- pba$model
	coefficients.star <- pba$coefficients.star
	coefficients.hat <- pba$coefficients.hat
	coefficients.hat.random <- pba$coefficients.hat.random
	
	# Create summary list to store results
	summary <- list()
	
	# Summarize estimate and confidience limits from observed model
	estimate <- coefficients.star[,'Estimate']
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
	
	# Apply transformation
	if (!is.null(transformation))
	{
		summary <- lapply(summary, function(x)
				{
					x[,-4] <- do.call(transformation, args=list(x[,-4], ...))
					return(x)
				})
	}
	
	# Calculate precision
	if (scale == "additive")
	{
		summary <- lapply(summary, function(x)
				{
					x$precision <- x[,3] - x[,2]
					return(x[ , c(1,2,3,5,4)])
				})
	}
	
	if (scale == "multiplicative")
	{
		summary <- lapply(summary, function(x)
				{
					x$precision <- x[,3] / x[,2]
					return(x[ , c(1,2,3,5,4)])
				})
	}
	
	# Return summary
	class(summary) <- "summary.pba"
	return(summary)
}




# 
pbaIterateSelection <- function(model, model.original=NULL, bias.tables, iter)
{
	# Create list to store updated models
	models.updated <- list()
	
	# Iterations
	for (i in 1:iter)
	{
		# Create copy of model to use during iterations
		model.updated <- model
		data <- model.updated$data
		rows <- nrow(data)
		
		selection <- llply(bias.tables, function(x)
				{
					x$selection
				})
		
		probs <- list()
		
		# Iterate variables
		for (j in names(selection))
		{	
			# Return if no selection bias defined for variable
			if (length(selection[[j]]) > 0)
			{			
				# Create empty vector to store selection probabilities
				prob <- vector(length = rows)
				
				# Find rows for a1s, a0s, b1s, and b0s
				rows.a1 <- which(data[,1]==1 & data[,j]==1)
				rows.a0 <- which(data[,1]==1 & data[,j]==0)
				rows.b1 <- which(data[,1]==0 & data[,j]==1)
				rows.b0 <- which(data[,1]==0 & data[,j]==0)
				
				# Assign probability of being selected
				prob[rows.a1] <- selection[[j]]$s.a1
				prob[rows.a0] <- selection[[j]]$s.a0
				prob[rows.b1] <- selection[[j]]$s.b1
				prob[rows.b0] <- selection[[j]]$s.b0
				
				# Store result
				probs[[j]] <- prob
			}
		}
		
		# Calculate joint probability of not being selected
		not.selected.prob <- 1 - Reduce("*", probs)
		
		# Over sample observations
		if (length(not.selected.prob) > 0)
		{
			not.selected <- rbinom(rows, 1, not.selected.prob)	
			model.updated$data <- rbind(data, data[which(not.selected==1),])
			models.updated[[i]] <- model.updated
		} else
		{
			models.updated[[i]] <- model
		}
	}
	
	return(models.updated)
}


#' Function to define distributions for bias parameters.
#' 
#' @param distr A character vector of length one naming a random generation 
#' ditribution.
#' @param args Additional arguments to be passed to the distr function. The 
#' argument \code{n} does not need to be specified, as it will be autumatically 
#' added when appropriate in subsequent functions.
#' 
#' @return An object of class \code{pba.distr}.
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

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




# Print method for pba objects
print.pba <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
	cat("\nA pba (Probabilistic Bias Analysis) object\n")
	cat(paste("\nIterations: ", x$iter, "\n\n", sep = ""))
	if (length(x$coefficients.star)) {
		cat("Original coefficients without bias adjustment")
		cat(":\n")
		print.default(format(x$coefficients.star, digits = digits), 
				print.gap = 2, quote = FALSE)
	}
	else cat("No coefficients\n\n")
	invisible(x)
}


# Print method for pba.variables object
print.pba.variables <- function(x, ...)
{
	if(length(x) == 1)
	{
		variables <- paste(": ", names(x), sep="")
	} else if(length(x) > 1)
	{
		variables <- paste("s: ", paste(names(x), collapse=", "), sep="")
	} else
	{
		cat("\n\nNo bias parameters defined]n")
		break()
	}
	cat(paste("\nA pba.variables object describing bias parameters for the following
							variable", variables, "\n", sep=""))
	for (i in names(x))
	{
		cat(paste("\nVariable ", i, ":\n", sep=""))
		if(!is.null(x[[i]]$misclassification[[1]]))
		{
			cat("\tMisclassification with:\n")
			cat(paste("\t\tSensitivity among cases distribution of", 
							x[[i]]$misclassification$se.a.distr$distr, "\n"))
			cat(paste("\t\tSensitivity among non-cases distribution of", 
							x[[i]]$misclassification$se.b.distr$distr, "\n"))
			cat(paste("\t\tSpecificity among cases distribution of", 
							x[[i]]$misclassification$sp.a.distr$distr, "\n"))
			cat(paste("\t\tSpecificity among non-cases distribution of", 
							x[[i]]$misclassification$sp.b.distr$distr, "\n"))		
		}
		if(!is.null(x[[i]]$selection[[1]]))
		{
			cat("\tSelection with\n")
			cat(paste("\t\tSelection among exposed cases distribution of", 
							x[[i]]$selection$s.a1.distr$distr, "\n"))
			cat(paste("\t\tSelection among non-exposed cases distribution", 
							x[[i]]$selection$s.a0.distr$distr, "\n"))
			cat(paste("\t\tSelection among exposed non-cases distribution of", 
							x[[i]]$selection$s.b1distr$distr, "\n"))
			cat(paste("\t\tSelection among non-exposed non-cases distribution of", 
							x[[i]]$selection$s.b0.distr$distr, "\n"))	
		}
		if(!is.null(x[[i]]$confounding[[1]][[1]]))
		{
			for (j in names(x[[i]]$confounding))
			{
				cat(paste("\tUnmeasured confounding by", x[[i]]$confounding[[j]]$name, "with:\n"))
				cat(paste("\t\tProportion among exposed distribution of",
								x[[i]]$confounding[[j]]$p1.distr$distr, "\n"))
				cat(paste("\t\tProportion among non-exposed distribution of",
								x[[i]]$confounding[[j]]$p0.distr$distr, "\n"))
				if(!is.null(x[[1]]$confounding[[j]]$rr.distr))
				{
					cat(cat(paste("\t\tRelative risk of confounding distribution of",
											x[[i]]$confounding[[j]]$rr.distr$distr, "\n")))
				}
				if(!is.null(x[[1]]$confounding[[j]]$rd.distr))
				{
					cat(cat(paste("\t\tRisk difference of confounding distribution of",
											x[[i]]$confounding[[j]]$rr.distr$distr, "\n")))
				}
			}
		}			
	}
}


#' Plot density ditribution of bias parameters of a pba object
#' 
#' @param pba A pba object.
#' @param density Logical; if true, plots smoothed density. If false, plots
#' histogram.
#' @param sclaes Character string passed to ggplot. This parameter controls
#' whether the x and y axes are free or fixed. Acceptable values are "free",
#' "fixed", "free_x", or "free_y"
#' @param types A character vector specifying the types of biases to plot. If
#' NULL, the function will plot all defined biases. Acceptable values are
#' "misclassification", "selection", and "confounding".
#' @param print Logical; if true, function prints the plots. If false, the
#' function only returns the plot objects in a list.
#' 
#' @return The plotBias function returns density plots or histograms of the
#' distribution of bias parameters for a pba object. If print = TRUE, the
#' function attempts to organize all plots into a single window. Individual 
#' plots may be accessed from the returned list. plotBias internally calls
#' pbaPlotMisclassification, pbaPlotSelection, pbaPlotConfoundingProportions,
#' and pbaPlotConfoundingRisks.
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# pba plot method
plotBias <- function(pba, density=T, scales='free', types=NULL, print=T)
{
	
	# Choose types of bias to plot
	if (is.null(types))
	{
		types <- laply(pba$bias.tables, function(x)
				{
					names(x)
				})
	}
	
	types <- unique(types)
	
	# Melt data into long-form
	data <- melt(pba$bias.tables)
	columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
	data <- data[ , columns]
	
	# Create plots list to store plots
	plots <- list()
	
	# Misclassification
	if ('misclassification' %in% types)
	{
		plots$misclassification <- pbaPlotMisclassification(data = data, 
				density = density, scales = scales)
	}
	
	# Selection
	if ('selection' %in% types)
	{
		plots$selection <- pbaPlotSelection(data = data, 
				density = density, scales = scales)
	}
	
	# Confounding proportions
	if ('confounding' %in% types)
	{
		plots$confounding.proportions <- pbaPlotConfoundingProportions(data = data, 
				density = density, scales = scales)
	}
	
	# Confounding risks
	if ('confounding' %in% types)
	{
		plots$confounding.risks <- pbaPlotConfoundingRisks(data = data, 
				density = density, scales = scales)
	}
	

	# Print plots with viewports
	if (print)
	{
		# Set up viewports
		n <- length(types)
		vp <- list()
		j <- 1
		for(i in types)
		{
			vp[[i]] <- viewport(x = (j - 0.5) / n, width = 1 / n, 
					just = c("centre", "centre"))
			j <- j + 1
		}
		
		# Print plots
		grid.newpage()
		pushViewport(viewport())
		popViewport()
		for(i in names(plots))
		{
			if (i %in% c("misclassification", "selection"))
			{
				pushViewport(vp[[i]])
				grid.draw(ggplotGrob(plots[[i]]))
				popViewport()
			} else if (i %in% c("confounding.proportions", "confounding.risks"))
			{
				pushViewport(vp[['confounding']])
				if (i == "confounding.proportions")
				{
					pushViewport(viewport(y = 1/3, height = 2/3, just=c("centre", "bottom")))
					grid.draw(ggplotGrob(plots[[i]]))
					popViewport()
				}
				if ( i == "confounding.risks")
				{
					pushViewport(viewport(y = 1/3, height = 1/3, just=c("centre", "top")))
					grid.draw(ggplotGrob(plots[[i]]))
					popViewport()
				}
				popViewport()
			}
		}		
	}
	
	# Return plots and viewports
	#result <- list(plots=plots, vps=vps)
	invisible(plots)
}


# Misclassification
pbaPlotMisclassification <- function(pba=NULL, data=NULL, density=T,
		parameters=c('se.a', 'se.b', 'sp.a', 'sp.b'), title='Misclassification', 
		scales='free', L2='misclassification')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
	
	# Subset misclassification data and parameters
	data <- subset(data, L2 == L2 & variable %in% parameters)
	
	# Define facet object
	facet <- facet_grid(variable~L1, scales = scales)
	
	# Define ggplot
	plot1 <- ggplot(data, aes(x=value)) + facet
	
	# Density
	plot1.density <- plot1 + geom_density() 
	
	# Histogram
	plot1.histogram <- plot1 + geom_histogram()
	
	if (density)
	{
		plot.result <- plot1.density
	} else
	{
		plot.result <- plot1.histogram
	}
	
	return(plot.result + opts(title = title))
}


# Confounding proportions
pbaPlotConfoundingProportions <- function(pba=NULL, data=NULL,
		parameters=c('p1', 'p0'), density=T, title='Confounding (proportions)', 
		scales='free', L2='confounding')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
	
	# Subset confounding data and proportions data
	data <- subset(data, L2 == L2 & !is.na(value) & 
					variable %in% parameters)
	
	# Define faet object
	facet <- facet_grid(variable~L3, scales=scales)
	
	## p1 and p0
	plot1 <- ggplot(data, aes(x=value)) + facet
	
	# Density
	plot1.density <- plot1 + geom_density() 
	
	# Histogram
	plot1.histogram <- plot1 + geom_histogram()
	
	if (density)
	{
		plot.result <- plot1.density
	} else
	{
		plot.result <- plot1.histogram
	}
	
	return(plot.result + opts(title = title))
}


# Confounding risks
pbaPlotConfoundingRisks <- function(pba=NULL, data=NULL,
		parameters=c('rr', 'rd'), density=T, title='Confounding (relative risks)',
		scales='free', L2='confounding')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
	
	# Subset confounding data and risk data
	data <- subset(data, L2 == L2 & !is.na(value) & 
					variable %in% parameters)
	
	# Define faet object
	facet <- facet_grid(variable~L3, scales=scales)
	
	## rr and rd
	plot1 <- ggplot(data, aes(x=value)) + facet
	
	# Density
	plot1.density <- plot1 + geom_density() 
	
	# Histogram
	plot1.histogram <- plot1 + geom_histogram()
	
	if (density)
	{
		plot.result <- plot1.density
	} else
	{
		plot.result <- plot1.histogram
	}
	
	return(plot.result + opts(title = title))
}



#' Plot density ditribution of estimates after adjusting for bias
#' 
#' @param pba A pba object.
#' @param data ?
#' @param exp Logical; if true, estimates are exponentiated. 
#' @param density Logical; if true, plots smoothed density. If false, plots
#' histogram.
#' @param sclaes Character string passed to ggplot. This parameter controls
#' whether the x and y axes are free or fixed. Acceptable values are "free",
#' "fixed", "free_x", or "free_y"
#' @param types A character vector specifying the types of biases to plot. If
#' NULL, the function will plot all defined biases. Acceptable values are
#' "misclassification", "selection", and "confounding".
#' @param print Logical; if true, function prints the plots. If false, the
#' function only returns the plot objects in a list.
#' 
#' @return The plotBias function returns density plots or histograms of the
#' distribution of bias parameters for a pba object. If print = TRUE, the
#' function attempts to organize all plots into a single window. Individual 
#' plots may be accessed from the returned list. plotBias internally calls
#' pbaPlotMisclassification, pbaPlotSelection, pbaPlotConfoundingProportions,
#' and pbaPlotConfoundingRisks.
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# Plot distribution of simulated estimates
plotEstimates <- function(pba=NULL, data=NULL, density=T, exp=F, adjust=1, 
		binwidth=NULL, scales='free', variables=NULL, print = T)
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
	
	# Print plot
	if (print)
	{
		print(plot)
	}
	
	# Return plot
	invisible(plot)
}



# Labeller
pbaLabeller <- function(x, y)
{
	if (y == 'se.a') c('Sensitivity', 'Sensitivity', 'Specificity', 'Specificity') 
	else y
}

pbaLabeller <- function(x, y)
{
	if (any(y %in% c('se.a', 'se.b', 'sp.a', 'sp.b')))
	{
		print(y)
		y[which(test=='se.a')] <- 'Sensitivity Outcome+'
		y[which(test=='se.b')] <- 'Sensitivity Outcome-'
		y[which(test=='sp.a')] <- 'Specificty Outcome+'
		y[which(test=='sp.b')] <- 'Specificity Outcome-'
	}
	
	return(y)
}



# Selection
pbaPlotSelection <- function(pba=NULL, data=NULL, density=T,
		parameters=c('s.a1', 's.a0', 's.b1', 's.b0'), title='Selection', 
		scales='free', L2='selection')
{
	pbaPlotMisclassification(pba = pba, data = data, density = density,
			parameters = parameters, title = title, scales = scales, L2 = L2)
}


#' Random generation from a trapezoidal distribution.
#' @param n Number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param min Lower bound of the trapezoid.
#' @param mode1 Lower mode of the trapezoid.
#' @param mode2 Upper mode of the trapezoid.
#' @param max Upper bound of the trapezoid.
#' 
#' @return A vector numbers randomly generated from the trapezoidal distribution.
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}
#' @references Adapted from Matthew Fox and colleagues' SAS macro SENSITIVITY ANALYSIS MISCLASSIFICATION MACRO version 1.1.
#' Available: \link{http://sites.google.com/site/biasanalysis/}

rtrapezoid <- function (n, min = 0, mode1 = 0.33, mode2 = 0.67, max = 1) 
{
	# Check arguments are valid
	if (length(n) > 1) 
		n <- length(n)
	if (n < 1 | is.na(n)) 
		stop(paste("invalid argument: n =", n))
	n <- floor(n)
	if (any(is.na(c(min, mode1, mode2, max)))) 
		return(rep(NaN, times = n))
	if (min > max | min > mode1 | min > mode2 | mode1 > max | mode2 > max |
			mode1 > mode2) 
		return(rep(NaN, times = n))
	if (any(is.infinite(c(min, mode1, mode2, max)))) 
		return(rep(NaN, times = n))
	
	# Generate random numbers
	p <- runif(n)
	r <- (p * (max + mode2 - min - mode1) + (min + mode1)) / 2
	
	# Adjust if random number is less than mode1
	rs.lt.mode1 <- which(r < mode1)
	r[rs.lt.mode1] <- min + sqrt((mode1 - min) * (2 * r[rs.lt.mode1] - min - mode1))
	
	# Adjust if random number is greater than mode2
	rs.gt.mode2 <- which(r > mode2)
	r[rs.gt.mode2] <- max - sqrt(2 * (max - mode2) * (r[rs.gt.mode2] - mode2))
	
	return(r)
}


#' Define bias parameters
#' 
#' @param variable A character vector naming the variable upon which the bias 
#' is effecting.
#' @param misclassification A list describing misclassification bias. The list 
#' consists of the following parameters:
#' \item{se.a.distr}{A pba.distr object defining the sensitivity among cases.}
#' \item{sp.a.distr}{A pba.distr object defining the specificity among cases.}
#' \item{se.b.distr}{A pba.distr object defining the sensitivity among 
#' non-cases.}
#' \item{sp.b.distr}{A pba.distr object defining the sensitivity among 
#' non-cases.}
#' \item{se.cor}{Correlation between sensitivity among cases and non-cases. 
#' Correlation of 1 indicates non-differential selection bias. Correlation of 
#' 0 indicates independent differential selection bias. Correlation less than 
#' 1 but greater than 0 indicates partial differential selection bias.}
#' \item{se.cor}{Correlation between specificity among cases and non-cases. 
#' Correlation of 1 indicates non-differential selection bias. Correlation of 0 
#' indicates independent differential selection bias. Correlation less than 1 
#' but greater than 0 indicates partial differential selection bias.}
#' @param selection A list describing selection bias. The list consists of the 
#' following parameters:
#' \item{s.a1.distr}{A pba.distr object defining selection among exposed cases.}
#' \item{s.a0.distr}{A pba.distr object defining selection among non-exposed 
#' cases.}
#' \item{s.b1.distr}{A pba.distr object defining selection among exposed 
#' non-cases.}
#' \item{s.b0.distr}{A pba.distr object defining selection among non-exposed 
#' non-cases.}
#' @param confounding A list containing one or more lists describing unmeasured 
#' confounding bias.  The confounder bias lists contain the following parameters:
#' \item{p1.distr}{A pba.distr object defining the probability of the 
#' unmeasured confounder among the exposed.}
#' \item{p0.distr}{A pba.distr object defining the probability of the 
#' unmeasured confounder among the non-exposed.}
#' \item{rr.distr}{(optional) A pba.distr object defining the relative risk 
#' association between the confounder and the outcome. rr.distr or or.distr must 
#' be defined. If both are defined, rr.distr is used instead of or.distr.}
#' \item{or.distr}{(optional) A pba.distr object defining the odds ratio 
#' association between the confounder and the outcome. rr.distr or or.distr 
#' must be defined. If both are defined, rr.distr is used instead of or.distr.}
#' 
#' @return The pba function returns an object of class \code{pba.variables}.
#' 
#' @export 
#' @author Jeremy Thoms Hetzel \email{jthetzel@@gmail.com}

# Function to create a define bias of a variable
pbaVariable <- function(variable, misclassification = NULL, 
		selection = NULL, confounding = NULL)
{											
	result <- list()
	result[[variable[[1]]]] <- list(variable=variable, 
			misclassification=misclassification, 
			selection=selection,
			confounding=confounding)
	
	class(result) <- "pba.variables"
	return(result)	
}






#' @aliases plotEstimates
#' @S3method plot pba

# Plot method for pba objects. Passes object to plotEstimates
plot.pba <- function(pba=NULL, data=NULL, density=T, exp=F, adjust=1, 
		binwidth=NULL, scales='free', variables=NULL, print = T)
{
	plotEstimates(pba = pba, data = data, density = density, exp = exp,
			adjust = adjust, binwidth = binwidth, scales = scales,
			variables = variables, print = print)
}
