# TODO: Add comment
# 
# Author: Антон
###############################################################################

p <- 1;

betta <- .4
gamma <- 10
delta <- .14



f <- function(x, z){
	c ( z[2], delta * z[1] * exp(gamma*betta*(1-z[1]) / (1 + betta * (1 - z[1])) ) )
	#c(z[2], 56*x^6)
}

K1 <- function(h, x.i, y.i){
	h*f(x.i, y.i)
}

K2 <- function(h, x.i, y.i){
	h*(f(x.i+h/2, y.i+K1(h, x.i, y.i)/2))
}

K3 <- function(h, x.i, y.i){
	h*f(x.i+h/2, y.i+K2(h, x.i, y.i)/2)		
}

K4 <- function(h, x.i, y.i){
	h*f(x.i+h, y.i+K3(h, x.i, y.i))
}

runge.kutta <- function(a, b, alpha ,n.c=100, print.graph=FALSE){
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	n <- length(x)
	y <- matrix(rep(0, times=2*n), n, 2)
	y[1,] <- c(alpha, 0)
	
	for (i in 1:(n-1))
		y[i+1,] <- it.r.k(h, x[i], y[i,])
	
	if (print.graph==TRUE){
		plot(x, y[,1], type="l", col="blue")
	}
	
	return (y[n, 1])
}

it.r.k <- function(h, x.i, y.i){
	y.i + 1/6 * (K1(h, x.i, y.i) + 2*K2(h, x.i, y.i) + 2*K3(h, x.i, y.i) + K4(h, x.i, y.i))
}

shooting.method <- function(epsilon=.001){
	a1 <- -10
	a2 <- -10
	it <- 100
	
	r.k.a1 <- runge.kutta(0, 1, a1) - 1
	r.k.a2 <- runge.kutta(0, 1, a2) - 1
	
	
	while (sign(r.k.a1) == sign(r.k.a2)){
		a2 <- a2 + .5
		r.k.a1 <- runge.kutta(0, 1, a1) - 1
		r.k.a2 <- runge.kutta(0, 1, a2) - 1
		it <- it - 1
		if (it == 0) stop("can't find alpha")
	}
	
	
	
	a3 <- (a2 + a1)/2
	r.k.a3 <- runge.kutta(0, 1, a3) - 1
	
	it <- 100
	while ( abs(r.k.a3) > epsilon ){
		
		if (sign(r.k.a1) != sign(r.k.a3)){
			a1 <- a1
			a2 <- a3
			r.k.a1 <- runge.kutta(0, 1, a1) - 1
			r.k.a2 <- runge.kutta(0, 1, a2) - 1
		}
		else{
			a1 <- a3
			a2 <- a2
			r.k.a1 <- runge.kutta(0, 1, a1) - 1
			r.k.a2 <- runge.kutta(0, 1, a2) - 1
		}
		
		a3 <- (a2 + a1)/2
		r.k.a3 <- runge.kutta(0, 1, a3) - 1
		
		it <- it - 1
		if (it == 0) stop("iteration process obviously won't converge")
	}
	
	return (a3)
}

shooting <- function(epsilon=.001){
	runge.kutta(0, 1, shooting.method(epsilon), print.graph=TRUE)
}

