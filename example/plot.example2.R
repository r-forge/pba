
exp <- pbaVariable(variable='exp',
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
								mode1=0.85, mode2=0.95, max=1.50)), 
				s.a0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b1.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				s.b0.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95))
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



exp2 <- pbaVariable(variable='exp2',
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
								mode1=0.85, mode2=0.95, max=1.50)), 
				s.a0.distr = pbaDistr('rtrapezoid', args=list(min=0.75, 
								mode1=0.85, mode2=0.95, max=1.00)),
				s.b1.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95)),
				s.b0.distr = pbaDistr('rtrapezoid', args=list(min=0.70, 
								mode1=0.80, mode2=0.90, max=0.95))
		),
		confounding = list(
				hypertension = list(
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


# Perform pba
iter <- 100
bias.tables <- pbaBiasTables(pba.variables = c(exp, exp2), iter = iter)
llply(bias.tables, function(x) { x$selection})
example$exp2 <- example$case
example$exp2[example$case==1] <- rbinom(example$exp2[example$case==1], 1, .55)
example$exp2[example$case==0] <- rbinom(example$exp2[example$case==0], 1, .45)
glm2 <- glm(case ~ exp + exp2, data=example, family=binomial())
pba6 <- pba(glm2, c(exp, exp2), iter=100)
summary(pba6)
pbaPlotEstimates(pba6)

pba6$bias.tables

# Try to figure out viewpoints
a <- pbaPlotBias(pba6, print=T)
b <- pbaPlotBias(pba6, print=T, types=c("misclassification", "selection", "confounding"))
c <- pbaPlotBias(pba6, print=T, types=c("selection", "confounding"))
d <- pbaPlotBias(pba6, print=T, types=c("misclassification", "confounding"))
e <- pbaPlotBias(pba6, print=T, types=c("misclassification", "selection"))
f <- pbaPlotBias(pba6, print=T, types=c("confounding", "selection", "misclassification", "selection"))

plot(pba6)
a <- plot(pba6, print = F)
plotBias(pba6)



vp1 <- viewport()
vp2 <- viewport(y = 1/3, height = 2/3, just=c("centre", "bottom"))
vp3 <- viewport(y = 1/3, height = 1/3, just=c("centre", "top"))

grid.newpage()
pushViewport(vp1)
popViewport()
pushViewport(vp2)
grid.draw(ggplotGrob(a$confounding.proportions))
popViewport()
pushViewport(vp3)
grid.draw(ggplotGrob(a$confounding.risks))
popViewport()

types <- c('misclassification', 'selection', 'confounding')

n <- length(types)

vp <- list()
j <- 1
for(i in types)
{
	vp[[i]] <- viewport(x = (j - 0.5) / n, width = 1 / n, 
			just = c("centre", "centre"))
	j <- j + 1
}


grid.newpage()
pushViewport(viewport())
popViewport()
for(i in names(a))
{
	if (i %in% c("misclassification", "selection"))
	{
		pushViewport(vp[[i]])
		grid.draw(ggplotGrob(a[[i]]))
		popViewport()
	} else
	{
		pushViewport(vp[['confounding']])
		if (i == "confounding.proportions")
		{
			pushViewport(viewport(y = 1/3, height = 2/3, just=c("centre", "bottom")))
			grid.draw(ggplotGrob(a[[i]]))
			popViewport()
		}
		if ( i == "confounding.risks")
		{
			pushViewport(viewport(y = 1/3, height = 1/3, just=c("centre", "top")))
			grid.draw(ggplotGrob(a[[i]]))
			popViewport()
		}
		popViewport()
	}
}
