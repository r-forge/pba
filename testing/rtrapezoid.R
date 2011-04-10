rtrapezoid <- function(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1)
{
	out <- qtrapezoid(p = runif(n), min = min, mode1 = mode1, mode2 = mode2,
			max = max, n1 = n1, n3 = n3, alpha = alpha)
	
	return(out)
}


alpha <- 1.2
n1 <- 1.5
n3 <- 2.5
r1 <- rtrapezoid(10000, n1=n1, n3=n3, alpha=alpha)
x <- seq(0,1,0.01)

plot(density(r1))
points(x,dtrapezoid(x, n1=n1, n3=n3, alpha=alpha))



alpha <- 5
ptrapezoid(qtrapezoid(.2, alpha=alpha), alpha=alpha)
