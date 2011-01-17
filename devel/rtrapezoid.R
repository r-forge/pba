# Function to randomly sample from a trapezoidal distribution
# Adapted from Fox M and colleagues'  
# SAS macro SENSITIVITY ANALYSIS MISCLASSIFICATION MACRO version 1.1.
# Available: 
rtrapezoid <- function (n = 1, min = 0, mode1 = 0.33, mode2 = 0.67, max = 1) 
{
	# Check arguments are valid
	if (length(n) > 1) 
		n <- length(n)
	if (n < 1 | is.na(n)) 
		stop(paste("invalid argument: n =", n))
	n <- floor(n)
	if (any(is.na(c(min, mode1, mode2, max)))) 
		return(rep(NaN, times = n))
	if (min > max | min > mode1 | min > mode2 | mode1 > max | mode2 > max |
			mode1 > mode2) 
		return(rep(NaN, times = n))
	if (any(is.infinite(c(min, mode1, mode2, max)))) 
		return(rep(NaN, times = n))
	
	# Generate random numbers
	p <- runif(n)
	r <- (p * (max + mode2 - min - mode1) + (min + mode1)) / 2
	
	# Adjust if random number is less than mode1
	rs.lt.mode1 <- which(r < mode1)
	r[rs.lt.mode1] <- min + sqrt((mode1 - min) * (2 * r[rs.lt.mode1] - min - mode1))
	
	# Adjust if random number is greater than mode2
	rs.gt.mode2 <- which(r > mode2)
	r[rs.gt.mode2] <- max - sqrt(2 * (max - mode2) * (r[rs.gt.mode2] - mode2))
	
	return(r)
}
