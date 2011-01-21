# Example of pba analysis
# Data from Fox et al
# 
# Author: jthetzel
###############################################################################

## Recreate Fox example analysis

# Require packages
require('plyr')
require('ggplot2')
require('QRMlib')

# Set working directory
directory <- "C:/Users/jthetzel/Research/pba"
setwd(directory)

# Load pba function
source('C:/Users/jthetzel/Research/pba/devel/pba.R')

# Set seed for reproducibility
set.seed(1234)

# Read in Fox's example set
example <- read.csv("C:/Users/jthetzel/Research/pba/other/SAS/example.sas7bdat.csv")

# Specify model
glm1 <- glm(case ~ exp, data=example, family=binomial())

# Specify bias variable, non-differential
exp.non.differential <- pbaVariable(variable='exp', 
												 se.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
												 							mode1=0.85, mode2=0.95, max=1.00)), 
												 sp.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00)),
												 se.b.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00)),
												 sp.b.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00)),
												 se.cor = 1,
												 sp.cor = 1)

# Specify bias variable, differential
exp.differential <- pbaVariable(variable='exp', 
												 se.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
																 mode1=0.85, mode2=0.95, max=1.00)), 
												 sp.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
																 mode1=0.85, mode2=0.95, max=1.00)),
												 se.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
																 mode1=0.80, mode2=0.90, max=0.95)),
												 sp.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
																 mode1=0.80, mode2=0.90, max=0.95)),
												 se.cor = 0.8,
												 sp.cor = 0.8)										 
							 
# Perform pba analysis, non-differential
pba1 <- pba(glm1, exp.non.differential, iter=1000, alpha=0.05)
lapply(pba1$summary, exp) # Original report in Fox:
													# Odds ratio (95% CL): 2.4 (1.2, 14)
pbaPlotBias(pba1)
pbaPlotBias(pba1, density=T)
pbaPlotEstimates(pba1)
pbaPlotEstimates(pba1, density=T)


# Perform pba analysis, differential
pba2 <- pba(glm1, exp.differential, iter=300, alpha=0.05)
lapply(pba2$summary, exp) # Original report in Fox:
													# Odds ratio (95% CL): 3.6 (1.6, 52)
pbaPlotBias(pba2)
pbaPlotEstimates(pba2)
pbaPlotEstimates(pba2, density=T)




# Create second exposure variable exp2
example$exp2 <- example$case
example$exp2[example$case==1] <- rbinom(example$exp2[example$case==1], 1, .55)
example$exp2[example$case==0] <- rbinom(example$exp2[example$case==0], 1, .45)

# Specify model
glm2 <- glm(case ~ exp + exp2, data=example, family=binomial())

# Specify second bias variable
exp2 <- pbaVariable(variable='exp2', 
		se.a.distr = pbaDistr('rnorm', args=list(mean=0.8, sd=0.05)),
		sp.a.distr = pbaDistr('rnorm', args=list(mean=0.9, sd=0.05)),
		se.b.distr = pbaDistr('rnorm', args=list(mean=0.8, sd=0.05)),
		sp.b.distr = pbaDistr('rnorm', args=list(mean=0.9, sd=0.05)),
		se.cor = 0.8,
		sp.cor = 0.9)


# Perform bias analysis with two bias variables
pba3 <- pba(glm2, c(exp.non.differential, exp2), iter=1000, alpha=0.05)
lapply(pba3$summary, exp)
pbaPlotBias(pba3)
pbaPlotBias(pba3, density=T)
pbaPlotEstimates(pba3)
pbaPlotEstimates(pba3, density=T)
pbaPlotEstimates(pba3, density=T, transform='exp')

# Specify poisson model
glm3 <- glm(case ~ exp + exp2, data=example, family=poisson())

# Perform bias analysis with two bias variables on poisson model
pba4 <- pba(glm3, c(exp.differential, exp2), iter=1000, alpha=0.05)
lapply(pba4$summary, exp)
pbaPlotBias(pba4)
pbaPlotBias(pba4, density=T)
pbaPlotEstimates(pba4)
pbaPlotEstimates(pba4, density=F)


















