# Fun with copulas
###############################################################################
require(mvtnorm)
require(lattice)

n <- 1000;
rho <- matrix(c(1, .4, .2, .4, 1, -.8, .2, -.8, 1), nrow = 3)

z <- rmvnorm(n, sigma = rho, method = "svd")
cor(z)

u <- pnorm(z, 0, 1)


##
##
n <- 5000

## create uncorrelated observations
se.a <- rtrapezoid(n, min=0.6, mode1=0.65, mode2=0.8, max=0.85)
se.b <- rtrapezoid(n, min=0.6, mode1=0.65, mode2=0.8, max=0.85)
sp.a <- rnorm(n, mean=0.9, sd=0.05)
sp.b <- rnorm(n, mean=0.9, sd=0.05)
#uncorrelated <- data.frame(se.a, se.b, sp.a, sp.b)
uncorrelated <- matrix(c(se.a, se.b, sp.a, sp.b), ncol=4)
row.names <- c('se.a','se.b','sp.a','sp.b')
colnames(uncorrelated) <- row.names

melted <- melt(uncorrelated)
densityplot(~ value | X2, data = melted, plot.points=F)



## correlation matrix
rho <- diag(4)
row.names <- c('se.a','se.b','sp.a','sp.b')
rownames(rho) <- colnames(rho) <- row.names
rho['se.a','se.b'] <- 0.9
rho['se.b','se.a'] <- 0.9
rho['sp.a','sp.b'] <- 0.7
rho['sp.b','sp.a'] <- 0.7

## create correlated observations from normal distributions
z <- rmvnorm(n, sigma = rho, method = "svd")

## De-normalize the correlated random draws with the normal moments 
## of the uncorrelated sensitivities and specificities
deNormalize <- function(normalized, original)
{
	sds <- apply(original, 2, sd)
	means <- apply(original, 2, mean)
	normalized * sds + means
}

correlated <- deNormalize( uncorrelated, z)

cor(correlated)

melted <- melt(correlated)
densityplot(~ value | X2, data = melted, plot.points=F)

#cor(pba2$bias.tables$exp$misclassification[,1:4], method="kendall")








## Copula
require(copula)
require(trapezoid)
require(ggplot2)

hists <- function(data)
{
	melted <- melt(data)
	if(is.null(melted$variable))
	{
		melted$variable <- melted$X2
	}
	out <- densityplot(~ value | variable, data = melted, plot.points=F)
	return(out)
}



myCop2 <- normalCopula(c(0, 0, 0), dim = 3, dispstr = "un")
myCop4 <- normalCopula(c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4), dim = 4, dispstr = "un")
myCop5 <- normalCopula(c(1, 0, 0, 0, 0, 1), dim = 4, dispstr = "un")
myCop6 <- normalCopula(c(0), dim = 2, dispstr = "un")


myMvd2 <- mvdc(copula = myCop2, margins = c("trapezoid", "trapezoid",
				"trapezoid"), paramMargins = list(list(0.7, 0.75, 0.85, 0.9), 
				list(0.7, 0.75, 0.85, 0.9), list(0.7, 0.75, 0.85, 0.9)))
myMvd4 <- mvdc(copula = myCop4, margins = c("trapezoid", "trapezoid",
				"trapezoid", "trapezoid"), paramMargins = list(list(0.6, 0.65, 0.75, 0.8), 
				list(0.7, 0.75, 0.85, 0.9), list(0.8, 0.85, 0.95, 1.0), list(0.8, 0.85, 0.95, 1.0)))
myMvd5 <- mvdc(copula = myCop5, margins = c("trapezoid", "trapezoid",
				"trapezoid", "trapezoid"), paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1)))
myMvd6 <- mvdc(copula = myCop6, margins = c("trapezoid", "trapezoid"), 
		paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1)))


a <- rmvdc(myMvd2, n = 10000)
hists(a)
cor(a)

b <- rmvdc(myMvd4, 10000)
hists(b)
cor(b)

random5 <- rmvdc(myMvd5, 10000)
hists(random5)
cor(random5)

random6 <- rmvdc(myMvd6, 10000)
hists(random6)
cor(random6)




myCop2 <- claytonCopula(c(.1, 0, 0, .9), dim = 4)
myCop2 <- frankCopula(c(.1, 0, 0, .9), dim = 4)
myCop2 <- huslerReissCopula(c(.1, 0, 0, .9))

myCop2 <- fgmCopula(c(.1, 0, 0, .9), dim=4)

myCop2 <- plackettCopula(c(0.5,0.5))

myMvd2 <- mvdc(copula = myCop2, margins = c("trapezoid", "trapezoid"), paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1)))
myMvd2 <- mvdc(copula = myCop2, margins = c("trapezoid", "trapezoid",
				"trapezoid", "trapezoid"), paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1)))

myCop2 <- normalCopula(c(.15643), dim = 2)
myMvd2 <- mvdc(copula = myCop2, margins = c("trapezoid", "trapezoid"), paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1)))

myCop2 <- normalCopula(c(-.1, .2, .3, .4, .5, .6), dim = 4, dispstr = "un")
myMvd2 <- mvdc(copula = myCop2, margins = c("trapezoid", "trapezoid",
				"trapezoid", "trapezoid"), paramMargins = list(list(0.75, 0.85, 0.95, 1), 
				list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1), list(0.75, 0.85, 0.95, 1)))

random1 <- rmvdc(myMvd2, 100000)
cor(random1)
cor(random1, method="kendall")
hists(random1)

random2 <- rmvdc(myMvd2, 1000)
cor(random2)
hists(random2)



fgm.cop <- fgmCopula(c(-.5,.5))
x <- rcopula(fgm.cop, 1000)
cor(x)


