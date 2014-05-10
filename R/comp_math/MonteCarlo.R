# TODO: Numerical Analysis. Work 10. Monte Carlo
# 
# Author: Антон
###############################################################################

eps <- .1

# no idea why it's does not work when we compute result using mean of 
# function vaues. and when we compute result by coefficient between 
# points we drop to shape area and points we drop out of shape area ???
getn <- function(x, a, b){		
	for(i in seq(from=30000, to=35000, by=100)){
		if (3*sqrt(pvar(x)/length(x)) > eps)
			x <- runif(i, a, b)
		else break
	}
	x
}

# variance for population (cause var(x) computes varience for sample 1/(n-1))
pvar <- function(population){
	mean(population^2) - mean(population)^2
}

rh <- .0001
hs <- .0001

rdev1 <- function(f, x0){
	(f(x0+rh/2) - f(x0-rh/2))/ rh
}

parint <- function(f1, f2, t0){
	f1(t0)*rdev1(f2, t0)
}

Simpson <- function(f, f1, f2, a, b, h=hs){
	x <- seq(from=a, to=b, by=h)
	n <- length(x)
	sum(abs(f(f1, f2, x[1:(n-1)])) + abs(f(f1, f2, x[2:n])) + 4*abs(f(f1, f2, x[2:n]-h/2)))*h/6
}

# try compute everything without epsilon precision
MC <- function(f, f1, f2, a, b){
	t <- runif(30000, a, b)
	y <- abs(f(f1, f2, t))
	(b-a)*mean(y)
}

splshape <- function(){
	t <- (1:4)
	x <- c(2, 5, 7, 2)
	y <- c(2, 7, 3, 2)
	a <- splinefun(t, x)
	b <- splinefun(t, y)
	plot(c(x,0), c(y,0), xlab="xaxis", ylab="yaxix")
	t <- c(seq(from=1, to=4, by=.025), 1)
	lines(a(t), b(t), col="red")
	print(Simpson(parint, a, b, 1, 3) - Simpson(parint, a, b, 3, 4))
	print(MC(parint, a, b, 1, 3) - MC(parint, a, b, 3, 4))
}

circle <- function(r=1){
	t <- seq(from=0, to=6.28, by=.02)
	plot(c(-3,3), c(-3, 3), xlab="xaxis", ylab="yaxix")
	lines(r*cos(t), r*sin(t), col="purple")
	print("Simpson")
	print(Simpson(parint, cos, sin, 0, 6.28))
	print("Monte Carlo #1")
	print(MC(parint, cos, sin, 0, 6.28))
	N <- 30000
	K <- 0
	x <- runif(N, -1, 1)
	y <- runif(N, -1, 1)
	for (i in 1:N){
		if (sqrt(x[i]^2+y[i]^2)<=r){
			K <- K+1
			putdot(x[i], y[i], "red")
		}
		else putdot(x[i], y[i], "blue")
	}
	print("Monte Carlo 2")
	print(2*2*K/N)
}

putdot <- function(x, y, color){
	for (i in 1:length(x))
		lines(c(x[i], x[i]+.001), c(y[i], y[i]+.001), col=color)
}
