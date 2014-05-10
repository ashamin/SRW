# TODO: Numeric Analysis. Work 8. Derivative computing
# 		Use after formDataDev.R
#		Interpolation function int is in functions.R
# 
# Author: Антон
###############################################################################


#error of function value
noise <- rnorm(length(x), mean=0, sd=0)

rh <- .01

rdev1 <- function(x0){
	(f(x0+rh/2) - f(x0-rh/2))/ rh
}

rdev2 <- function(x0){
	(f(x0+rh) - 2*f(x0) + f(x0-rh)) / (rh^2)
}

rdev3 <- function(x0){
	(f(x0+3*rh/2) - 3*f(x0+rh/2) + 3*f(x0-rh/2) - f(x0-3*rh/2)) / (rh^3)
}

dev1 <- function(x0){
	n <- length(y)
	res  <- (y[3:n]-y[1:(n-2)])/(2*h)
	res + (noise[3:n]-noise[1:(n-2)])/(2*h)
}

dev2 <- function(x0){
	n <- length(y)
	res <- (y[5:n]-2*y[3:(n-2)]+y[1:(n-4)])/(4*h^2)
	res + (noise[5:n]-2*noise[3:(n-2)]+noise[1:(n-4)])/(4*h^2)
}

dev3 <- function(x0){
	n <- length(y)
	res <- (y[7:n]-3*y[5:(n-2)]+3*y[3:(n-4)]-y[1:(n-6)])/(8*h^3)
	res + (noise[7:n]-3*noise[5:(n-2)]+3*noise[3:(n-4)]-noise[1:(n-6)])/(8*h^3)
}


dv1 <- function(){
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1 <- matrix(rep(0, times=4*length(criteria(x, x))), length(criteria(x, x)), 4)
	d1[,1] <- criteria(x, x) 	
	d1[,2] <- criteria(rdev1(x), x)
	d1[,3] <- criteria(c(NA, dev1(x), NA), x)
	#error of function value
	noise <<- gnoise
	d1[,4] <- criteria(c(NA, dev1(x), NA), x)
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1
}	


dv2 <- function(){
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1 <- matrix(rep(0, times=4*length(criteria(x, x))), length(criteria(x, x)), 4)
	d1[,1] <- criteria(x, x) 	
	d1[,2] <- criteria(rdev2(x), x)
	d1[,3] <- criteria(c(NA, NA, dev2(x), NA, NA), x)
	#error of function value
	noise <<- gnoise
	d1[,4] <- criteria(c(NA, NA, dev2(x), NA, NA), x)
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1
}	

dv3 <- function(){
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1 <- matrix(rep(0, times=4*length(criteria(x, x))), length(criteria(x, x)), 4)
	d1[,1] <- criteria(x, x) 	
	d1[,2] <- criteria(rdev3(x), x)
	d1[,3] <- criteria(c(NA, NA, NA, dev3(x), NA, NA, NA), x)
	#error of function value
	noise <<- gnoise
	d1[,4] <- criteria(c(NA, NA, NA, dev3(x), NA, NA, NA), x)
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1
}

#old version sample
initdv3 <- function(){
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1 <- matrix(rep(0, times=4*length(x)), length(x), 4)
	d1[,1] <- x 	
	d1[,2] <- c(rdev3(x))
	d1[,3] <- c(NA, NA, NA, dev3(x), NA, NA, NA)
	#error of function value
	noise <<- gnoise
	d1[,4] <- c(NA, NA, NA, dev3(x), NA, NA, NA)
	noise <<- rnorm(length(x), mean=0, sd=0)
	d1
}	

