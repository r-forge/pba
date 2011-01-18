# TODO: Add comment
# 
# Author: jthetzel
###############################################################################

# Set seed for reproducibility
set.seed(1234)

# Read in Fox's example set
example <- read.csv("C:/Users/jthetzel/Research/sensitivity/SAS/example.sas7bdat.csv")

# Create second exposure variable
example$exp2 <- sample(c(0,1), nrow(example), replace=T)

# Logistic regression with single exposure
glm1 <- glm(case ~ exp, data=example, family=binomial())

# Logistic regression with two exposures
glm2 <- glm(case ~ exp + exp2, data=example, family=binomial())

# Create pbaDistr using triangle distributions
require(triangle)
se.a.distr.triangle <- pbaDistr(distr='rtriangle', args=list(a=.8, b=.9, c=.85))
sp.a.distr.triangle <- pbaDistr(distr='rtriangle', args=list(a=.8, b=.9, c=.85))
se.b.distr.triangle <- pbaDistr(distr='rtriangle', args=list(a=.8, b=.9, c=.85))
sp.b.distr.triangle <- pbaDistr(distr='rtriangle', args=list(a=.8, b=.9, c=.85))

# Create pbaDistr and using normal distributions
se.a.distr.norm <- pbaDistr(distr='rnorm', args=list(0.7, 0.05))
sp.a.distr.norm <- pbaDistr(distr='rnorm', args=list(0.8, 0.05))
se.b.distr.norm <- pbaDistr(distr='rnorm', args=list(0.9, 0.05))
sp.b.distr.norm <- pbaDistr(distr='rnorm', args=list(0.7, 0.1))

# Create two pbaVariables
exp <- pbaVariable(variable='exp', se.a.distr.triangle, sp.a.distr.triangle, 
									 se.b.distr.triangle, sp.b.distr.triangle,
									 se.cor=.8, sp.cor=.9)
exp2 <- pbaVariable(variable='exp2', se.a.distr.norm, sp.a.distr.norm, 
									 se.b.distr.norm, sp.b.distr.norm,
									 se.cor=.9, sp.cor=.8)

# Test pba function with single bias variable
time0 <- proc.time()
test <- pba(glm1, exp, iter=100)
time1 <- proc.time()
time <- time1 - time0

# Test pba with multiple bias variables
biases <- c(exp, exp2)
time0 <- proc.time()
test <- pba(glm2, biases, iter=100)
time1 <- proc.time()
time <- time1 - time0



# Try fixing unreal values
replaced <- list()

iter <- 1000
i <- 'exp'
model <- glm1
table.star <- table(model$model[,1], model$model[,i])
a1.star <- table.star[4]
a0.star <- table.star[2]
b1.star <- table.star[3]
b0.star <- table.star[1]
a1.stars <- rep(a1.star, iter)
a0.stars <- rep(a0.star, iter)
b1.stars <- rep(b1.star, iter)
b0.stars <- rep(b0.star, iter)

correlated.tables <- pbaBiasCor(exp, iter)
current.table <- correlated.tables[[i]]

result <- pbaBackCalculate(a1.star=a1.stars, 
														   a0.star=a0.stars, 
															 b1.star=b1.stars,
															 b0.star=b0.stars,
															 se.a = current.table$se.as,
															 sp.a = current.table$sp.as,
															 se.b = current.table$se.bs,
															 sp.b = current.table$sp.bs)

#result$or.hat <- with(result, (a1.hat / a0.hat) / (b1.hat / b0.hat))	

apply(result,2,min)
apply(result,2,max)	
which(result < 0)


rows <- which(apply(result[,c('a1.hat', 'a0.hat', 'b1.hat', 'b0.hat')], 1, function(x)
{
	any(x < 0)
}))
result[rows,]

correlated.table.replacement <- pbaBiasCor(exp, length(rows))

dim(correlated.tables[[i]][rows,])
dim(correlated.table.replacement[[i]])

correlated.tables[[i]][rows,] <- correlated.table.replacement[[i]]

result.replacement <- pbaBackCalculate(a1.star=a1.stars[rows], 
		a0.star=a0.stars[rows], 
		b1.star=b1.stars[rows],
		b0.star=b0.stars[rows],
		se.a = correlated.tables[[i]]$se.as[rows],
		sp.a = correlated.tables[[i]]$sp.as[rows],
		se.b = correlated.tables[[i]]$se.bs[rows],
		sp.b = correlated.tables[[i]]$sp.bs[rows])

result[rows,] <- result.replacement

if (is.null(replaced[[i]]))
{
	replaced[[i]] <- 0
}
replaced[[i]] <- replaced[[i]] + length(rows)




# Try to model probability distribution?
r <- rtrapezoid(1000, .75, .8, .9, .95)
plot(density(r))




length(pba1$coefficients.hat$exp[,'Estimate'][pba1$coefficients.hat$exp[,'Estimate']<0])/
length(pba1$coefficients.hat$exp[,'Estimate'])

exp(min(pba1$coefficients.hat$exp[,'Estimate']))




test <-  matrix(c(seq(1:81)), nrow=9, byrow=T)
test2 <- matrix(c(rep(9,27)), nrow=3, byrow=T)

test[c(1,6,9),] <- test2[,]

temp <- pbaBiasCor(exp.differential, iter=1+1)
temp <- pbaMisclassification(glm1, exp.differential, iter=1)

pba(glm1, exp.differential, iter=10, alpha=0.05)


hist(pba1$coefficients.hat.random$exp)

# Histograms
variables <- names(pba1$coefficients.hat.random)
variables <- 'exp'
data <- melt(pba1$coefficients.hat.random)
data <- subset(data, L1 %in% variables)

lower <- quantile
p1 <- ggplot(data, aes(x=value))
p2 <- p1 + geom_histogram(binwidth=0.1, drop=T) + xlim(quantile(value, .01), quantile(value, .99))
p2
p2 <- p1 + geom_density(trim=T) + facet_grid(.~L1, scales='free', space='free')
p3 <- p2
p3
#histogram
data <- melt(pba1$coefficients.hat.random[['exp']])
p1 <- ggplot(data, aes(x=value))
p2 <- p1 + geom_histogram(binwidth=0.1)
p3 <- p2 + xlim(quantile(data$value, 0.001), quantile(data$value, 0.999))
p3
#density
data <- melt(pba1$coefficients.hat.random[['exp']])
p1 <- ggplot(data, aes(x=value))
p2 <- p1 + geom_density()
p3 <- p2 + xlim(quantile(data$value, 0.001), quantile(data$value, 0.999))
p3 


# Biases
data <- melt(pba1$bias$tables$exp[c('se.a', 'se.b', 'sp.a', 'sp.b')])
data <- melt(pba2$bias$tables$exp[c('se.a', 'se.b', 'sp.a', 'sp.b')])
data <- melt(pba2$bias$tables)

p1 <- ggplot(data, aes(x=value))
p2 <- p1 + geom_density() + facet_grid(variable~.)
p2


temp <- melt(pba2$bias$tables)

temp <- lapply(pba2$bias$tables, function(x)
		{
			x[c('se.a', 'se.b', 'sp.a', 'sp.b')]
		})


















while (length(rows >= 1))
{
	# Create replacement sensitivities and specificities
	correlated.table.replacement <- pbaBiasCor(pba.variables[i], length(rows))
	
	# Replace original sensitivities and specificities with new values
	correlated.tables[[i]][rows,] <- correlated.table.replacement[[i]]
	
	# Back calculate expected counts from new sensitivities and specificities
	result.replacement <- pbaBackCalculate(a1.star=a1.stars[rows], 
			a0.star=a0.stars[rows], 
			b1.star=b1.stars[rows],
			b0.star=b0.stars[rows],
			se.a = correlated.tables[[i]]$se.as[rows],
			sp.a = correlated.tables[[i]]$sp.as[rows],
			se.b = correlated.tables[[i]]$se.bs[rows],
			sp.b = correlated.tables[[i]]$sp.bs[rows])
	
	# Replace rows in result table with new values
	print('begin')
	print(result[,])
	result[rows,] <- result.replacement
	print(result[,])
	
	# Record number of rows replaced
	if (is.null(replaced[[i]]))
	{
		replaced[[i]] <- 0
	}
	replaced[[i]] <- replaced[[i]] + length(rows)
	
	# Identify rows with counts less than 0
	rows <- which(apply(result[,c('a1.hat', 'a0.hat', 'b1.hat', 'b0.hat')], 1, 
					function(x)
					{
						any(x < 0)
					}))
	print(rows)
}