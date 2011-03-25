
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












a <- pbaPlotBias(pba6, print=F)

n <- length(a)
vps <- vpList()
for (i in 1:n)
{
	vps[[i]] <- viewport(x=(i-.5)/n, y=1/2, w=1/n, h=1, name=names(a)[i], just=c('center'))
}

vp <- vpTree(parent=viewport(w=1, h=1),
				children=vps
		)


grid.newpage()
pushViewport(vp)		
for(i in names(a))
{
	seekViewport(i)
	grid.draw(ggplotGrob(a[[i]]))
}


# Grob object for the two confounding objects
conf <- gTree(ggplotGrob(a$confounding.proportions), ggplotGrob(a$confounding.risks), name='confounding')

vp.conf.p <- viewport(x=1/2, y=2/3, w=1, h=2/3)
vp.conf.r <- viewport(x=1/2, y=1/3, w=1, h=1/3)
conf.p <- grob(ggplotGrob(a$confounding.proportions), vp=vp.conf.p)
conf.r <- grob(ggplotGrob(a$confounding.risks), vp=vp.conf.r, name='risks')
conf.p <- ggplotGrob(a$confounding.proportions)
conf.r <- ggplotGrob(a$confounding.risks)

grid.newpage()
pushViewport(vp)
grid.draw(conf.p)
grid.draw(conf.r)

conf.tree <- grobTree(conf.p, conf.r)


