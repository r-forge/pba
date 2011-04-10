# TODO: Add comment
# 
# Author: jthetzel
###############################################################################

sum.exp <- function(x)
{
	x <- summary(x)
	out <- as.data.frame(exp(x$coefficients[,-4]))
	out$p <- x$coefficients[,4]
	return(out)
}



## Specify distributions of variables
n <- 10000
variables <- list()
variables$crc <- pbaDistr("binom", size = 1, prob = 0.05)
variables$male <- pbaDistr("binom", size = 1, prob = 0.5) #http://www.google.com/url?sa=t&source=web&cd=1&ved=0CBQQFjAA&url=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FList_of_countries_by_sex_ratio&ei=16qgTczKH4StgQe01e3aBQ&usg=AFQjCNEqUhvWXHHduD3LdvcmJTubA-yCww&sig2=kwqlzPTzcbJlliEvIC4yoQ
variables$smoke <- pbaDistr("binom", size = 1, prob = 0.208) #http://www.cdc.gov/mmwr/preview/mmwrhtml/mm5644a2.htm
variables$colonoscopy <- pbaDistr("binom", size = 1, prob = 0.45)
variables$history <- pbaDistr("binom", size = 1, prob = 0.15)

correlation <- diag(1, length(variables))
rownames(correlation) <- colnames(correlation) <- names(variables)

correlation['crc', 'male'] <- 0.1
correlation['crc', 'smoke'] <- 0.3
correlation['crc', 'colonoscopy'] <- -0.4
correlation['crc', 'history'] <- 0.3
correlation['male', 'smoke'] <- 0.2
correlation['male', 'colonoscopy'] <- -0.2
correlation['male', 'history'] <- 0.0
correlation['smoke', 'colonoscopy'] <- -0.1
correlation['smoke', 'history'] <- 0.05
correlation['colonoscopy', 'history'] <- 0.2

correlation[lower.tri(correlation)] <- t(correlation)[lower.tri(t(correlation))]
sigma <- correlation[lower.tri(correlation)]

data <- pbaCorrelate(n, variables, sigma)
cor(data)

data <- data.frame(data)

glm1 <- glm(crc~male+smoke+colonoscopy+history, data=data, family=poisson())
sum1 <- summary(glm1)
sum1.exp <- sum.exp(glm1)




## Induce bias
data.biased <- data
# history is unknown confounder
data.biased <- data.biased[-which(names(data.biased) == "history")]
# smoking misclassified
# se.a=


