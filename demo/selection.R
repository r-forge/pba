exp.sel <- pbaVariable(variable='exp',
		misclassification = list( 
				se.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.50)), 
				sp.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				se.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				sp.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				se.cor = 0.8,
				sp.cor = 0.8
		),
		selection = list(
				s.a1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.a0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00))			
		)
)

exp2.sel <- pbaVariable(variable='exp2',
		misclassification = list( 
				se.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.50)), 
				sp.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				se.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				sp.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				se.cor = 0.8,
				sp.cor = 0.8
		),
		selection = list(
				s.a1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.a0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00))			
		)
)


# Variables
iter <- 1

# Calculate bias tables from pba.variables
bias.tables <- pbaBiasTables(pba.variables = c(exp.sel, exp2.sel), iter = iter)

#
selection <- llply(bias.tables, function(x)
		{
			x$selection
		})

data <- glm2$model
data[which(data$case==1 & data$exp==1),]$exp * (1 - selection$exp$s.a1)

probs <- list()

prob <- vector(length=nrow(data))
prob[data$case==1 & data$exp==1] <- selection$exp$s.a1
prob[data$case==1 & data$exp==0] <- selection$exp$s.a0
prob[data$case==0 & data$exp==1] <- selection$exp$s.b1
prob[data$case==0 & data$exp==0] <- selection$exp$s.b0
probs[['exp']] <- prob

prob <- vector(length=nrow(data))
prob[data$case==1 & data$exp2==1] <- selection$exp2$s.a1
prob[data$case==1 & data$exp2==0] <- selection$exp2$s.a0
prob[data$case==0 & data$exp2==1] <- selection$exp2$s.b1
prob[data$case==0 & data$exp2==0] <- selection$exp2$s.b0
probs[['exp2']] <- prob

not.selected.probs <- 1 - Reduce("*", probs)

not.selected <- rbinom(not.selected.probs, 1, not.selected.probs)

resample <- data[which(not.selected==1),]
data.new <- rbind(data, resample)

selected <- vector(length = nrow(data.new))
selected[1:nrow(data)] <- T


test <- pbaIterateSelection(model = glm2, model.original=NULL, bias.tables = bias.tables, iter = 1)
	
test2 <- update(test[[1]], data=test[[1]]$data)


# Test
pba5 <- pba(glm2, c(exp.sel, exp2.sel), iter=100, alpha=0.05)
lapply(pbaSummary(pba5), exp) # Original report in Fox:
# Odds ratio (95% CL): 2.4 (1.2, 14)
pbaPlotEstimates(pba5)
pbaPlotSelection(pba=pba5)
pbaPlotBias(pba = pba5, types = 'selection')

