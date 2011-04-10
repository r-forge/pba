# Example of pba analysis
# Data from Fox et al
# 
# Author: jthetzel
###############################################################################

## Recreate Fox example analysis

# Require packages
require('pba')
require('copula')
require('ggplot2')
require('trapezoid')
require('fUtilities')

# Set working directory
directory <- "C:/Users/jthetzel/Research/pba"
setwd(directory)
directory <- "/home/jthetzel/Research/pba/"
setwd(directory)

# Load pba function
#source('C:/Users/jthetzel/Research/pba/pkg/R/pba.R')
#source('C:/Users/jthetzel/Research/pba/pkg/R/plots.R')
#source('C:/Users/jthetzel/Research/pba/pkg/R/rtrapezoid.R')
#source('/home/jthetzel/Research/pba/pkg/R/pba.R')
#source('/home/jthetzel/Research/pba/pkg/R/plots.R')
#source('/home/jthetzel/Research/pba/pkg/R/rtrapezoid.R')

# Set seed for reproducibility
set.seed(1234)

# Read in Fox's example set
example <- read.csv("C:/Users/jthetzel/Research/pba/other/SAS/example.sas7bdat.csv")
example <- read.csv("/home/jthetzel/Research/pba/other/SAS/example.sas7bdat.csv")

# Specify model
glm1 <- glm(case ~ exp, data=example, family=binomial())

# Specify bias variable, non-differential
exp.non.differential <- pbaVariable(variable='exp',
												 misclassification = list(
												 se.a.distr = pbaDistr('trapezoid', min=0.75, 
												 							mode1=0.85, mode2=0.95, max=1.00), 
												 sp.a.distr = pbaDistr('trapezoid', min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00),
												 se.b.distr = pbaDistr('trapezoid', min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00),
												 sp.b.distr = pbaDistr('trapezoid', min=0.75, 
																			mode1=0.85, mode2=0.95, max=1.00),
												 se.cor = 1,
												 sp.cor = 1))

# Specify bias variable, differential
exp.differential <- pbaVariable(variable='exp',
	misclassification = list( 
		se.a.distr = pbaDistr('trapezoid', min=0.75, 
			mode1=0.85, mode2=0.95, max=1.00), 
	 	sp.a.distr = pbaDistr('trapezoid', min=0.75, 
			mode1=0.85, mode2=0.95, max=1.00),
	 	se.b.distr = pbaDistr('trapezoid', min=0.70, 
			mode1=0.80, mode2=0.90, max=0.95),
	 	sp.b.distr = pbaDistr('trapezoid', min=0.70, 
			mode1=0.80, mode2=0.90, max=0.95),
	 	se.cor = 0.8,
	 	sp.cor = 0.8,
	 	sigma = c(0.8, 0, 0, 0, 0, 0.8)))										 
								 
# Perform pba analysis, non-differential
pba1 <- pba(glm1, exp.non.differential, iter=100)
summary(pba1, transformation="exp", scale="multiplicative") 
													# Original report in Fox:
													# Odds ratio (95% CL): 2.4 (1.2, 14)
plotBias(pba1)
pbaCor(pba1)
plotEstimates(pba1, variables=c(-1))
plotEstimates(pba1, variables=c(-1), exp=T)
plotEstimates(pba1, density=F)
a <- plotEstimates(pba1, density=F, exp=T)
a + xlim(0,5)
plotEstimates(pba1, density=F, binwidth=0.1)
plotEstimates(pba1, exp=T)


# Perform pba analysis, differential
pba2 <- pba(glm1, exp.differential, iter=1000)
summary(pba2, transformation="exp", scale="multiplicative")  # Original report in Fox:
													# Odds ratio (95% CL): 3.6 (1.6, 52)
plotBias(pba2)
plotBias(pba2, density=F)
plotEstimates(pba2)
plotEstimates(pba2, adjust=0.8)
plotEstimates(pba2, exp=T)
plotEstimates(pba2, exp=T, adjust=0.6)
plotEstimates(pba2, density=F)




# Create second exposure variable exp2
example$exp2 <- example$case
example$exp2[example$case==1] <- rbinom(example$exp2[example$case==1], 1, .55)
example$exp2[example$case==0] <- rbinom(example$exp2[example$case==0], 1, .45)

# Specify model
glm2 <- glm(case ~ exp + exp2, data=example, family=binomial())

# Specify second bias variable
exp2 <- pbaVariable(variable='exp2', 
		misclassification = list(
			se.a.distr = pbaDistr('norm', args=list(mean=0.8, sd=0.05)),
			sp.a.distr = pbaDistr('norm', args=list(mean=0.9, sd=0.05)),
			se.b.distr = pbaDistr('norm', args=list(mean=0.8, sd=0.05)),
			sp.b.distr = pbaDistr('norm', args=list(mean=0.9, sd=0.05)),
			se.cor = 0.8,
			sp.cor = 0.9)
)


# Perform bias analysis with two bias variables
pba3 <- pba(glm2, c(exp.non.differential, exp2), iter=100)
summary(pba3, transformation="exp", scale="multiplicative") 
plotBias(pba3)
plotBias(pba3, density=F)
plotEstimates(pba3)
plotEstimates(pba3, density=F)
plotEstimates(pba3, exp=T)

# Specify poisson model
glm3 <- glm(case ~ exp + exp2, data=example, family=poisson())

# Perform bias analysis with two bias variables on poisson model
pba4 <- pba(glm3, c(exp.differential, exp2), iter=100)
summary(pba4, transformation="exp", scale="multiplicative")
plotBias(pba4)
plotBias(pba4, density=T)
plotEstimates(pba4)
plotEstimates(pba4, density=F)
plotEstimates(pba4, exp=T)
plotEstimates(pba4, exp=T, variables=c(-1))


# Adjust for misclassification and selection bias
exp.ms <- pbaVariable(variable='exp', 
		misclassification = list(
				se.a.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.a.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.b.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.b.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.cor = 0.8,
				sp.cor = 0.9),
		selection = list(
				s.a1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.a0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				s.b1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.b0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				s.a.cor = 1, s.b.cor = 1, sigma=c(1,0,0,0,0,1))
)
exp2.ms <- pbaVariable(variable='exp2', 
		misclassification = list(
				se.a.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.a.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.b.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.b.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.cor = 0.8,
				sp.cor = 0.9),
		selection = list(
				s.a1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.a0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				s.b1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.b0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				sigma = c(0, 0, 0, 0, 0, 0))
)

pba5 <- pba(glm2, c(exp.ms, exp2.ms), iter=100)
summary(pba5, scale="multiplicative", transformation="exp")
plotEstimates(pba5)
plotBias(pba5, print=T)
plotBias(pba5, print=T, types="misclassification")
pbaCor(pba5)






# Adjust for misclassification, selection bias, and unmeasured confounding
exp.msc <- pbaVariable(variable='exp', 
		misclassification = list(
				se.a.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.a.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.b.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.b.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.cor = 0.8,
				sp.cor = 0.9),
		selection = list(
				s.a1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.a0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				s.b1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.b0.distr = pbaDistr('norm', mean=0.9, sd=0.05)),
		confounding = list(
				smoking = list(
						p1.distr = pbaDistr('trapezoid', min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00),
						p0.distr = pbaDistr('trapezoid', min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00),
						rr.distr = pbaDistr('trapezoid', min=2.5, 
										mode1=3.0, mode2=3.5, max=4.0),
						name = 'smoking'
				))
)
exp2.msc <- pbaVariable(variable='exp2', 
		misclassification = list(
				se.a.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.a.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.b.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				sp.b.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				se.cor = 0.8,
				sp.cor = 0.9),
		selection = list(
				s.a1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.a0.distr = pbaDistr('norm', mean=0.9, sd=0.05),
				s.b1.distr = pbaDistr('norm', mean=0.8, sd=0.05),
				s.b0.distr = pbaDistr('norm', mean=0.9, sd=0.05)),
		confounding = list(
			htn = list(
					p1.distr = pbaDistr('trapezoid', min=0.75, 
									mode1=0.85, mode2=0.95, max=1.00),
					p0.distr = pbaDistr('trapezoid', min=0.75, 
									mode1=0.85, mode2=0.95, max=1.00),
					rr.distr = pbaDistr('trapezoid', min=2.5, 
									mode1=3.0, mode2=3.5, max=4.0),
					p.cor = 0.8,
					name = 'hypertension'))
)


pba6 <- pba(glm2, c(exp.msc, exp2.msc), iter=1000)
summary(pba6, scale="multiplicative", transformation="exp")
plot(pba6)
plot(pba6, exp=T)
plotBias(pba6)







