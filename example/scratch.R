
# Specify bias variable, differential
exp.test <- pbaVariable(variable='exp',
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
		confounding = list(
				smoking = list(
						p1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00)),
						p0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00)),
						rr.distr = pbaDistr('rtrapezoid', args=list(min=2.5, 
										mode1=3.0, mode2=3.5, max=4.0)),
						name = 'smoking'
				)
		)
)


# Specify bias variable, differential
exp.test2 <- pbaVariable(variable='exp2',
		misclassification = list( 
				se.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.25)), 
				sp.a.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				se.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				sp.b.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				se.cor = 0.8,
				sp.cor = 0.8
		),
		confounding = list(
				htn = list(
						p1.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00)),
						p0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
										mode1=0.85, mode2=0.95, max=1.00)),
						rr.distr = pbaDistr('rtrapezoid', args=list(min=2.5, 
										mode1=3.0, mode2=3.5, max=4.0)),
						name = 'hypertension'
				)
		)
)


pba3 <- pba2(model=glm2, pba.variables=c(exp.test, exp.test2), iter=100, alpha=0.05)
pbaPlotEstimates(pba3)
pbaPlotBiasMisclassification(pba3)
pbaPlotBiasConfounding(pba3)
pbaSummary(pba3)




### plyr gglopt2
bias.tables <- pba3$bias.tables
parameters <- c('se.a', 'se.b', 'sp.a', 'sp.b', 'p1', 'p0', 'rr', 'rd')
a <- melt(bias.tables)
columns <- which(names(a) %in% c('variable', 'value', 'L1', 'L2', 'L3'))
a <- a[ , columns]
a <- subset(a, variable %in% parameters)
a <- subset(a, !is.na(value))

b <- subset(a, L2 == 'confounding')


head(subset(a, ))

p1 <- pbaPlotBiasMisclassification(pba3); p1
p1 <- pbaPlotBiasConfounding(pba3); p1

pbaPlotEstimates(pba3)