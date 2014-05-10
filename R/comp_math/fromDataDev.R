# TODO: Add comment
# 
# Author: Антон
###############################################################################


f <- function(x0){
	x0^4
}

#step
inith <- 0.1
h <- inith

#vector of arguments for table
x <- seq(from=0-6*h, to=0+56*h, by=h)
y <- f(x)

gnoise <- rnorm(length(x), mean=0, sd=0.2)

criteria <- function(sample, factor){
	#why don't work??
	sample[abs(as.integer(factor) - factor) < .00001]
	#sample[seq(from=7, to=(length(factor)-6), by=1/h)]
}
