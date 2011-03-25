# TODO: Add comment
# 
# Author: jthetzel
###############################################################################


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
pba6 <- pba(glm2, c(exp, exp2), iter=100, alpha=0.05)
pbaSummary(pba6)
pbaPlotEstimates(pba6)

pba6$bias.tables


pbaPlotBias(pba6)
plots <- list()
plots$misclassification <-pbaPlotMisclassification(pba6)
plots$selection <-pbaPlotSelection(pba6)
plots$confounding.proportions <- pbaPlotConfoundingProportions(pba6)
plots$confounding.risks <- pbaPlotConfoundingRisks(pba6)

g.pm <- ggplotGrob(pm)
grid.draw(g.pm)

vp <- viewport(width = 0.5, height = 0.2, name='vp')
pushViewport(vp)
grid.draw(g.pm)
popViewport()

n <- length(plots)
for(i in 1:n)
{
	vp <- viewport(width = 1/n, height = 1, x = (i/n) - (1/(2*n)))
	pushViewport(vp)
	grid.draw(ggplotGrob(plots[[i]]))
	popViewport()
}

g1 <- ggplotGrob(plots$confounding.proportions)
g2 <- ggplotGrob(plots$confounding.risks)

vp1 <- viewport(width = 1, height = 2/3, x = 1/2, y = 4/6)
vp2 <- viewport(width = 1, height = 1/3, x = 1/2, y = 1/6)
vp3 <- viewport(width = 1, height = 1, x = 1/2, y = 2/3)
vp4 <- viewport(width = 1, height = 1/3, x = 1/2, y = 1/3)

grid.newpage()
pushViewport(vp1)
grid.draw(g1)
popViewport()
pushViewport(vp2)
grid.draw(g2)
popViewport()

t1 <- gTree(g1, vp=vp1)
t2 <- gTree(g2, vp=vp2)
list1 <- gList(g2)
tt <- gTree(g1, vp=vp1, children=list1)

grid.draw(tt)



grid.newpage()
tree <- vpTree(viewport(w=1, h=1, name="Figures"),
		vpList(
				viewport(w=1/3, x=1/6, h=1, y=0.5, name="Misclassification"),
				viewport(w=1/3, x=1/2, h=1, y=0.5, name="Selection"),
				vpTree(viewport(w = 1/3, x = 5/6, h = 1, y = 0.5, name="Confounding"),
						vpList(
								viewport(x=1/2, y=2/3, w=1, h=2/3, 
										name="Proportions"),
								viewport(x=1/2, y=1/6, w=1, h=1/3, 
										name="Risk")))
				))
pushViewport(tree)
for (i in c("Figures", "Misclassification", "Selection", "Confounding", 
		"Proportions", "Risk")) {
	seekViewport(i)
	grid.rect()
	grid.text(current.vpTree(FALSE),
			x=unit(1, "mm"), y=unit(1, "npc") - unit(1, "mm"),
			just=c("left", "top"),
			gp=gpar(fontsize=8))
}

list1 <- glist()

g.p <- ggplotGrob(plots$confounding.proportions) 