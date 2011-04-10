
from sympy import *

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
n1 = Symbol('n1')
n3 = Symbol('n3')
alpha = Symbol('alpha')
p = Symbol('p')
q = Symbol('q')
x = Symbol('x')

one = solve((2 * alpha * n1 * n3) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1)) * ((x - a) / (b - a))**(n1 - 1) - p, x)

two = solve(((2 * n1 * n3) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * 
(c - b) * n1 * n3) + (2 * (d - c) * n1))) *
((alpha - 1) * ((c - x) / (c - b)) + 1) - p, x)

solve((2 * n1 * n3) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * 
(c - b) * n1 * n3) + (2 * (d - c) * n1)) *
((d - x) / (d - c))**(n3 -1) - p, x)

###
lower = nsolve((2 * alpha * (b - a) * n3) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1)) * ((q - a) / (b - a))**(n1) - p, q)
one = (2 * alpha * (b - a) * n3) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1)) * ((q - a) / (b - a))**(n1) - p

middle = solve(((2 * alpha * (b - a) * n3) + (2 * (q - b) * n1 * n3 * (1 + ((alpha - 1) / 2) * ((2*c - b - q) / (c - b))))) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1)) - p, q)

upper = solve(1 - ((2 * (d - c) * n1) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1))) * ((d - q) / (d - c))**(n3) - p, q)
three = (1 - ((2 * (d - c) * n1) / ((2 * alpha * (b - a) * n3) + ((alpha + 1) * (c - b) * n1 * n3) + (2 * (d - c) * n1))) * (((d - q) / (d - c))**(n3)) - p)
threeSimp = hypersimp(three, n3)


a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=0.5;
a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=0.1;
a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=0.9;
a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=1.5;
a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=3;
a=0;b=1/3;c=2/3;d=1;n1=2;n3=2;alpha=1;pi=3;

test <- function(pi, a=0, b=1/3, c=2/3, d=1, n1=2, n3=2, alpha=1)
{
lower <- (a**2*(b - a) + b**2*(b - a) - a**2*(b*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1)) - a*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1))) + a*b*(b*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1)) - a*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1))) - 2*a*b*(b - a))/((b - a)*(b*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1)) - a*(c*pi/2 - b*pi/2 + b*pi/n1 + c*pi/(2*alpha) - a*pi/n1 - b*pi/(2*alpha) + d*pi/(alpha*n3) - c*pi/(alpha*n3))**(1/(1 - n1))))

# two:
middle <- (-2*b*n1*n3 - 2*b*c*n1*pi - 2*c*d*n1*pi + 2*alpha*c*n1*n3 + 2*b*d*n1*pi - 2*a*alpha*b*n3*pi - 2*alpha*b*c*n3*pi + 2*a*alpha*c*n3*pi + 2*b*c*n1*n3*pi + 2*alpha*b*c*n1*n3*pi + 2*n1*pi*c**2 - n1*n3*pi*b**2 - n1*n3*pi*c**2 + 2*alpha*n3*pi*b**2 - alpha*n1*n3*pi*b**2 - alpha*n1*n3*pi*c**2)/(-2*n1*n3 + 2*alpha*n1*n3)

# three:
upper <- (d**2*(d*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3)) - c*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3))) - c**2*(d - c) - d**2*(d - c) - c*d*(d*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3)) - c*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3))) + 2*c*d*(d - c))/((d - c)*(d*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3)) - c*(c*pi/2 - b*pi/2 + alpha*c*pi/2 - alpha*b*pi/2 + d*pi/n3 - c*pi/n3 + alpha*b*pi/n1 - a*alpha*pi/n1)**(1/(1 - n3))))

return(cbind(lower, middle, upper))

}







## YACAS
OldSolve((2*alpha*(b-a)*n)/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m))*((q-a)/(b-a))^(m)==p,q)
OldSolve(((2*alpha*(b-a)*n)+(2*(q-b)*m*n*(1+((alpha-1)/2)*((2*c-b-q)/(c-b)))))/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m))==p,q)
OldSolve(1-((2*(d-c)*m)/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m)))*((d-q)/(d-c))^(n)==p,q)
# Try new Solve
Solve(((2*alpha*(b-a)*n)+(2*(q-b)*m*n*(1+((alpha-1)/2)*((2*c-b-q)/(c-b)))))/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m))==p,q)


#lower
a=0; b=1/3; c=2/3; d=1; m=2; n=2; alpha=1; q=0.1; p=0.0225
((2*alpha*(b-a)*n)/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m))*((q-a)/(b-a))^(m))
((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*alpha*(b-a)*n))^(1/m)*(b-a)+a

#middle
alpha=1;q=0.5;p=0.375
((2*alpha*(b-a)*n)+(2*(q-b)*m*n*(1+((alpha-1)/2)*((2*c-b-q)/(c-b)))))/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m))
2*c-b-(((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)/(2*(q-b)*m*n)-1)*2*(c-b))/(alpha-1)

Solve(2*c-b-(((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)/(2*(q-b)*m*n)-1)*2*(c-b))/(alpha-1), q)
b+((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)*2*(c-b))/(2*(2*c-b-((-2)*(c-b))/(alpha-1))*(alpha-1)*n*m)

#upper
alpha=1.5;q=0.75;p=0.8875
1-((2*(d-c)*m)/((2*alpha*(b-a)*n)+((alpha+1)*(c-b)*m*n)+(2*(d-c)*m)))*((d-q)/(d-c))^(n)
d-(((1-p)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*(d-c)*m))^(1/n)*(d-c)



























       