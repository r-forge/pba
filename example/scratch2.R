head(example)

# Impute values with given RR
rr <- 4
p.a <- .5
p.b <- .4

# rr <- (a1 / (a1 + b1)) / (b1 / (b1 + b0))
# p.a <- a1 / (a1 + a0)
# p.b <- b1 / (b1 + b0)

a1=105; a0=85; b1=527; b0=93
n <- a1 + a0 + b1 + b0

test <- pbaBackCalculateConfounding(a1=105, a0=85, b1=527, b0=93, 
		p1=0.80, p0=0.05, rr=0.63, rd=NULL)

test <- pbaConfounding(a1=105, a0=85, b1=527, b0=93, 
		p1=0.80, p0=0.05, rr=0.63, rd=NULL)

test.rd <- pbaConfounding(a1=105, a0=85, b1=527, b0=93, 
		p1=0.80, p0=0.05, rr=NULL, rd= -0.37)

p1 <- 0.80
p0 <- 0.05
rd <- -.037
t11 <- (a1 + b1) * p1
t01 <- (a0 + b0) * p0
t10 <- (a1 + b1) - t11
t00 <- (a0 + b0) - t10
t1  <- a1 + b1
t0  <- a0 + b0

a11 <- ((rd * t1 * (t1 - t11)) + (a1 * t11)) / t1
a01 <- ((rd * t0 * (t0 - t01)) + (a0 * t01)) / t0

RD <- -0.37
m <- 632
M1 <- 505.6
a <- 105

A1 <- (RD * m * (m - M1) + (a * M1)) / m


















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
			), 
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


test <- pbaSampleConfounding(exp.test$exp$confounding, 100)
names(test)
plot(density(test$smoking$p1))
max(test$smoking$p1)

pbaSampleMisclassification(exp.test$exp$misclassification, 100)

tables <- pbaBiasTables(exp.test, 100)

tables2 <- pbaCalculatePredictiveValues(tables, glm1, 100)

tables3 <- pbaCorrectMisclassification(tables2, exp.test, glm1, 100)






data <- glm1$model
data$smoking <- rbinom(nrow(data), 1, .5)
glm9 <- update(glm1, data = data, formula = '. ~ . + smoking')






#
pbaIterateMisclassification <- function(exposure, data, misclassification, 
		iter)
{
	
	# Find rows for a1s, a0s, b1s, and b0s
	rows.a1 <- which(data[,1]==1 & data[,exposure]==1)
	rows.a0 <- which(data[,1]==1 & data[,exposure]==0)
	rows.b1 <- which(data[,1]==0 & data[,exposure]==1)
	rows.b0 <- which(data[,1]==0 & data[,exposure]==0)
	
	# Calculate a1., a0., b1., and b0.
	a1. <- length(rows.a1)
	a0. <- length(rows.a0)
	b1. <- length(rows.b1)
	b0. <- length(rows.b0)
	n1. <- a1. + b1.
	n0. <- a0. + b0.
	
	# Calculate n11, n01, n10, and n00
	n11 <- n1. * p1
	n01 <- n0. * p0
	n10 <- n1. - n11
	n00 <- n0. - n01
	
	# Calculate a11 and a01
	a11 <- (rr * n11 * a1.) / ((rr * n11) + n1. - n11)
	a01 <- (rr * n01 * a0.) / ((rr * n01) + n0. - n01)
	
	
	# Simulate occurence of correct classification for a1s, a0s, b1s, and b0s
	correct.a1 <- rbinom(length(rows.a1), 1, misclassification$ppv.a[i])
	correct.a0 <- rbinom(length(rows.a0), 1, misclassification$npv.a[i])
	correct.b1 <- rbinom(length(rows.b1), 1, misclassification$ppv.b[i])
	correct.b0 <- rbinom(length(rows.b0), 1, misclassification$npv.b[i])
	
	# Change exposure if classification not correct
	data[rows.a1,exposure][correct.a1==0] <- 
			as.numeric(!data[rows.a1,exposure][correct.a1==0])
	data[rows.a0,exposure][correct.a0==0] <- 
			as.numeric(!data[rows.a0,exposure][correct.a0==0])
	data[rows.b1,exposure][correct.b1==0] <- 
			as.numeric(!data[rows.b1,exposure][correct.b1==0])
	data[rows.b0,exposure][correct.b0==0] <- 
			as.numeric(!data[rows.b0,exposure][correct.b0==0])
	
	return(data)
}









#########
#########
#########

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
				), 
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

model <- glm1
pbaVariables <- exp.test
iter <- 1

# Sample bias parameters
bias.tables <- pbaBiasTables(pbaVariables, iter = iter)
bias.tables$exp$misclassification
bias.tables

# Calculate predictive values for misclassification biases
bias.tables <- pbaCalculatePredictiveValues(bias.tables = bias.tables,
		model = model, iter = iter)
bias.tables
bias.tables$exp$misclassification

# Correct unrealistic misclassification values
bias.tables <- pbaCorrectMisclassification(bias.tables = bias.tables,
		pbaVariables = pbaVariables, model = model, iter = iter)
bias.tables$exp$misclassification

i <- pbaVariables$exp

misclassification <- pbaSampleMisclassification(
		misclassification = i$misclassification, iter = 1)

# Asdjust data for misclassification bias
data.adjusted <- pbaIterateMisclassification(model = model, 
		bias.tables = bias.tables, iter = 1)[[1]]

# Confounding
exposure <- 'exp'
outcome <- 'case'
model <- glm1
name <- 'htn'
confounding <- bias.tables$exp$confounding$hypertension
model.updated <- pbaIterateConfoundingSingle(exposure = exposure, 
		outcome = outcome, model = model,	confounding = confounding, name = name, 
		iter = 1)[[1]]
model.updated$formula


models.updated <- pbaIterateConfounding(model = model, 
		bias.tables = bias.tables, iter = 1)[[1]]
models.updated$formula
head(model.updated$data)

# pba2
test <- pba2(model = glm2, pba.variables = exp.test, iter = 100)






bt <- pbaBiasTables(c(exp.test, exp.test2), iter=1)




a <- pba3$bias.lists
bias.lists <- a
b <- pbaBiasListsTables(pba3$bias.lists)


confounding <- ldply(bias.lists[['exp']], function(x)
		{
			ldply(x$confounding, function(y)
					{
						y
					})
		})

confounding <- split(confounding, f = confounding$.id)





pba1$model.summaries[[2]]



a <- pbaPlotBiasConfounding(pba1)
b <- pbaPlotBiasConfounding(pba3)

vp1 <- viewport(width = 1, height = 2/3, x = 1/2, y = 4/6)
vp2 <- viewport(width = 1, height = 1/3, x = 1/2, y = 1/6)
print(b[[1]], vp = vp1)
print(b[[2]], vp = vp2)

p1 <- pbaPlotBias(pba3)
pbaPlotEstimates(pba3)
print(p1$plots[[1]], vp = p1$vps[[1]])
print(p1$plots[[2]], vp = p1$vps[[2]])
print(p1$plots[[3]], vp = p1$vps[[3]])







# check memory
test <- list()
for(i in 1:100000)
{
	test[[i]] <- glm1$data
}




Rprof()
pba1 <- pba(glm1, exp.non.differential, iter=100, alpha=0.05)
Rprof(NULL)
summaryRprof()
