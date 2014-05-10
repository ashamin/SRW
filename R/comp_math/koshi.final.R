# TODO: Add comment
# 
# Author: Антон
###############################################################################

p <- 1;

alpha <- c(0, 1/2, 1/2, 1)
betta <- matrix(c(rep (0, times=12)), 4, 3)
betta[2, ] <- c(1/2, 0, 0)
betta[3, ] <- c(0, 1/2, 0)
betta[4, ] <- c(0, 0, 1)
p.y <<- c(1/6, 1/3, 1/3, 1/6)
p.y.precision <- c(rep(0, times=4))

g <- function(x){
	x^3
	#1
}

gd <- function(x){
	3*x^2
}

f <- function(x, u){
	gd(x) + p*(u - g(x))
	#-30*u
}

k <- function(i, h, x.i, y.i){
	if (i != 1)
		for (j in (i-1):1)
			yy <- 
	else 
}

compute.n.plus.1.value <- function(h, x.i, y.i, q){
	res <- y.i
	for (i in 1:q)
		res <- res + p.y[i] * k(i, h, x.i, y.i)
	return (res)
}

runge.kutta <- function(a, b, n.c, auto.step=FALSE){
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	coefs <- classic.coef.init()
	n <- length(x)
	y <- c(rep(0, times=n))
	y[1] = g(x[1])
	for (i in 1:(n-1))
		y[i+1] <- compute.n.plus.1.value(h, x[i], y[i], 4)
	
	#par(mfcol = c(1, 2))
	plot(x, y, type="l", col="blue")
	title("Рунге-Кутта")
	#lines(x, g(x), type="l", col="blue")
	for (i in seq(from=1, to=length(x), by=floor(length(x)/10)))
		mprint(x[i], y[i])
}


Felberg.coef.init <- function(){
	alpha <<- c(0, 1/4, 3/8, 12/13, 1, 1/2)
	betta <<- matrix(c(rep (0, times=30)), 6, 5)
	betta[2, ] <<- c(1/4, rep(0, times=4))
	betta[3, ] <<- c(3/32, 9/32, rep(0, times=3))
	betta[4, ] <<- c(1932/2197, -7200/3197, 7296/2197, rep(0, times=2))
	betta[5, ] <<- c(439/216, -8, 3680/513, -845/4104, 0)
	betta[6, ] <<- c(-8/27, 2, -3544/2565, 1859/4104, -11/40)
	p.y <<- c(25/216, 0, 1408/2565, 2197/4104, -1/5, 0)
	p.y.precision <<- c(16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)
	list(a=alpha, b=betta, p.y=p.y, p.y.precision)
}

classic.coef.init <- function(){
	alpha <<- c(0, 1/2, 1/2, 1)
	betta <<- matrix(c(rep (0, times=12)), 4, 3)
	betta[2, ] <<- c(1/2, 0, 0)
	betta[3, ] <<- c(0, 1/2, 0)
	betta[4, ] <<- c(0, 0, 1)
	p.y <<- c(1/6, 1/3, 1/3, 1/6)
	p.y.precision <- c(rep(0, times=4))
	list(a=alpha, b=betta, p.y=p.y)
}

runge.kutta.adams.diff <- function(a, b, n.c){
	par(mfcol = c(1, 2))
	runge.kutta(a, b, n.c)
	adams(a, b, n.c)
}