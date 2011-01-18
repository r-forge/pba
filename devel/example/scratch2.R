# TODO: Add comment
# 
# Author: jthetzel
###############################################################################


# example of how to use a Guassian Copula to create random draws from
# normally distributed data

require("reshape")
require("plyr")
require("QRMlib")
require("Matrix")

#make this reproducable
set.seed(2)

#how many draws in our starting data set?
n <- 1e4

# how many draws do we want from this distribution?
drawCount <- 1e4

# ourData will be my starting data which we'll
# base our model on
myData <- rnorm(n, 5, .6)
yourData <- myData + rnorm(n, 8, .25)
hisData <- myData + rnorm(n, 6, .4)
herData <- hisData + rnorm(n, 8, .35)

myData <- rnorm(n, 5, .6)
yourData <- myData + rnorm(n, 8, .25)
hisData <- myData + rnorm(n, 6, .4)
herData <- hisData + rnorm(n, 8, .35)

ourData <- data.frame(myData, yourData, hisData, herData)

# now we have raw correlations in the 70s, 80s, and 90s. Funny how that
# works out
cor(ourData)

#set up some simple functions for normalizing (Z Score) and
#denormalizing
normalMoments <- function(t) {
	as.list(c(mean=mean(t), sd=sd(t)))
}

normalize <- function(x) {
	(x - mean(x)) / sd(x)
}

deNormalize <- function(x, sampleMean, sampleSd) {
	(x * sampleSd + sampleMean)
}


# create an object with the mean and sd of each
# of the margins... you can do this without plyr
# using apply() if you're into that sort of thing
# but alply makes it so easy to work with the results
# we will use these later to denormalize
normalMomentList <- alply(ourData, 2, normalMoments)

# normalize i.e. create Z score
# I think you don't HAVE to do this with gaussian marginals,
# but you DO with non-gaussian... so I'm in the habit.
ourDataZ <- data.frame(apply(ourData, 2, normalize))

# now ourDataZ is a data frame of Z scores.
# prove this to yourself by checking the
# mean and SD of each margin

apply(ourDataZ, 2, mean)
apply(ourDataZ, 2, sd)

#looks like mean of 0 and sd of 1 to me

# this fixes positive def, fits a copula, and makes draws
# in one line of code. Try that in Excel, bitches.
myDraws <- rcopula.gauss(n=drawCount,
		Sigma=as.matrix(nearPD(Spearman(ourDataZ),corr=TRUE)$mat),
		d=ncol(ourData))

# rcopula.gauss spits out points from 0-1 (i.e. q values) so we need to turn those into
# Z scores by doing a little norm inversion.
myDraws <- qnorm(myDraws)

#check the mean and sd
apply(myDraws, 2, mean)
apply(myDraws, 2, sd)

#should be 0, 1... and they are close..

# but we can do better than that... we know statistics!
# let's do a Kolmogorov-Smirnov test to see if these
# samples come from the same distribution
# KS does not check correlation,
# it only tests if two sets of samples
# came from same dist.. we'll check each column

for (i in 1:ncol(ourData)){
	print(ks.test(myDraws[[i]], ourDataZ[[i]]))
}

# if the p-value of the KS test is < .05 then we
# reject that the distributions are equal
# they all look > .05 to me.
# Kolmogorov-Smirnov makes me want to drink, for some odd reason
# If your starting sample is small, you'll notice that a couple of the variables
# fail or just barely pass the KS test. This is common because KS is non-parametric and
# the starting sample will be 'lumpy' and not that big. If starting n gets up
# to, say, 10K, then they all do better. That's why I started with 10K and still
# it can be hit or miss

# so we have a metric shit ton of random draws. But they are Z scores
# and we want them put back in their original shapes. So let's do that:

myDrawsDenorm <- myDraws
for (i in 1:ncol(myDrawsDenorm)) {
	myDrawsDenorm[,i] <- deNormalize(myDraws[,i],
			normalMomentList[[i]][[1]],
			normalMomentList[[i]][[2]])
}

myDrawsDenorm <- data.frame(myDrawsDenorm)
names(myDrawsDenorm) <- names(ourData)

#let's look at the mean and standard dev of the starting data
apply(ourData, 2, mean)
apply(ourData, 2, sd)

#compare that with our sample data
apply(myDrawsDenorm, 2, mean)
apply(myDrawsDenorm, 2, sd)


# so myDrawsDenorm contains the final draws
# let's check Kolmogorov-Smirnov between the starting data
# and the final draws

for (i in 1:ncol(ourData)){
	print(ks.test(myDrawsDenorm[[i]], ourData[[i]]))
}

# holy shit. it works!

#look at the correlation matrices
cor(myDrawsDenorm)
cor(ourData)

#it's fun to plot the variables and see if the PDFs line up
#you could do this for each variable. It's a good sanity check.
plot(density(myDrawsDenorm$myData))
lines(density(ourData$myData), col="red")

# there's a test to see if the corr matricices are the same
# but I'm too lazy to google for it


