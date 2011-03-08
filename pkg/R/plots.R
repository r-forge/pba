# pba plot method
pbaPlotBias <- function(pba, density=T, scales='free', print=T)
{
	# Melt data into long-form
	data <- melt(pba$bias.tables)
	columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
	data <- data[ , columns]
	
	# Misclassification
	plot.misclassification <- pbaPlotMisclassification(data = data, 
			density = density, scales = scales)
	
	# Confounding proportions
	plot.confounding.proportions <- pbaPlotConfoundingProportions(data = data, 
			density = density, scales = scales)
	
	# Confounding risks
	plot.confounding.risks <- pbaPlotConfoundingRisks(data = data, 
			density = density, scales = scales)
	
	# Define viewports
	vp.misclassification <- viewport(width = 1/2, height = 1, x = 1/4, y = 1/2)
	vp.confounding.proportions <- viewport(width = 1/2, height = 2/3, x =3/4, y = 4/6)
	vp.confounding.risks <- viewport(width = 1/2, height = 1/3, x =3/4, y = 1/6)
	
	# Print plots with viewports
	if (print)
	{
		print(plot.misclassification, vp = vp.misclassification)
		print(plot.confounding.proportions, vp = vp.confounding.proportions)
		print(plot.confounding.risks, vp = vp.confounding.risks)
	}
	
	# Save plots with viewports
	plots <- list(misclassification = plot.misclassification, 
			confounding.proportions = plot.confounding.proportions,
			confounding.risks = plot.confounding.risks)
	vps <- list(misclassification = vp.misclassification, 
			confounding.proportions = vp.confounding.proportions,
			confounding.risks = vp.confounding.risks)
	
	# Return plots and viewports
	result <- list(plots=plots, vps=vps)
	return(result)
}


# Misclassification
pbaPlotMisclassification <- function(pba=NULL, data=NULL, density=T,
		parameters=c('se.a', 'se.b', 'sp.a', 'sp.b'), title='Misclassification', 
		scales='free')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
	
	# Subset misclassification data and parameters
	data <- subset(data, L2 == 'misclassification' & variable %in% parameters)
	
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
		scales='free')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
	
	# Subset confounding data and proportions data
	data <- subset(data, L2 == 'confounding' & !is.na(value) & 
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
		scales='free')
{
	if (is.null(data))
	{
		# Melt data into long-form
		data <- melt(pba$bias.tables)
		columns <- which(names(data) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
		data <- data[ , columns]
	}
		
	# Subset confounding data and risk data
	data <- subset(data, L2 == 'confounding' & !is.na(value) & 
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


# Plot distribution of simulated estimates
pbaPlotEstimates <- function(pba=NULL, data=NULL, density=T, exp=F, adjust=1, 
		binwidth=NULL, scales='free', variables=NULL)
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