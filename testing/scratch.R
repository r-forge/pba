# TODO: Add comment
# 
# Author: jthetzel
###############################################################################

test <- function(p, a=0, b=1/3, c=2/3, d=1, m=2, n=2, alpha=1)
{
	
	one <- b+((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)*2*(c-b))/(2*(2*c-b-((-2)*(c-b))/(alpha-1))*(alpha-1)*n*m)
	two <- ((sqrt((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))^2/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)^2-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)))-(((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha))
	three <- (-((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)+sqrt((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))^2/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)^2-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha))

	one.true <- ptrapezoid(q = one, min=a, mode1=b, mode2=c, max=d, n1=m, n3=n, alpha=alpha)
	two.true <- ptrapezoid(q = two, min=a, mode1=b, mode2=c, max=d, n1=m, n3=n, alpha=alpha)
	three.true <- ptrapezoid(q = three, min=a, mode1=b, mode2=c, max=d, n1=m, n3=n, alpha=alpha)
	
	cbind(p, one, one.true, two, two.true, three, three.true)
}

p <- seq(0,1,.1)
test(p=p, alpha=2)


