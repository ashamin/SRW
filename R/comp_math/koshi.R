# TODO: Add comment
# 
# Author: ?????
###############################################################################

p <- 1;

g <- function(x){
	x^8
	#x^3
	#1
}

gd <- function(x){
	8*x^7
	#3*x^2
}

f <- function(x, u){
	gd(x) + p*(u - g(x))
	#-30*u
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

E <- function(h, x.i, y.i){
	2/3 * (K1(h, x.i, y.i) - K2(h, x.i, y.i) - K3(h, x.i, y.i) + K4(h, x.i, y.i))
}

runge.kutta <- function(a, b, n.c=100, auto.step=TRUE, epsilon=.00001){
	h <- (b-a)/n.c
	h.max <- (b-a)/5
	h.min <- (b-a)/20000
	x <- c(a)
	n <- length(x)
	y <- c(rep(0, times=n))
	y[1] = g(x[1])
	i <- 1
	
	if (auto.step==TRUE){
		err <- abs(E(h, x[i], y[i])) / max(abs(y[i]), 1)
		if (err > epsilon) #h <- h * min(1.5, max(.7, .9*(epsilon/err)^.2))
			h <- h/2
	}
	
	while (x[i] < b){
		if (auto.step == TRUE){
			err <- abs(E(h, x[i], y[i])) / max(abs(y[i]), 1)
			if (err < epsilon){
				
				y[i+1] <- it.r.k(h, x[i], y[i])
				x <- c(x, x[i] + h)
				i <- i + 1
			}
			else{
				h <- h/2
				#h <- h * min(1.5, max(.7, .9*(epsilon/err)^.2))
				#h <- h * max(h.min, min(h.max, .8*h*(epsilon/err)^.2))
			}
		}
		else{
			y[i+1] <- it.r.k(h, x[i], y[i])
			x <- c(x, x[i] + h)
			i <- i + 1
		} 
			
	}
	
	print(length(x))
	
	x <- x[1:(length(x)-1)]
	y <- y[1:(length(y)-1)]

	#for (i in 1:(n-1)){
	#	y[i+1] <- it.r.k(h, x[i], y[i])
	#}
	
	#par(mfcol = c(1, 2))
	plot(x, y, type="l", col="blue")
	title("Runge-Kutta")
	#lines(x, g(x), type="l", col="blue")
	for (i in seq(from=1, to=length(x), by=floor(length(x)/10)))
		mprint(x[i], y[i])
}

it.r.k <- function(h, x.i, y.i){
	y.i + 1/6 * (K1(h, x.i, y.i) + 2*K2(h, x.i, y.i) + 2*K3(h, x.i, y.i) + K4(h, x.i, y.i))
}

adams <- function(a, b, n.c){
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	n <- length(x)
	y <- c(rep(0, times=n))
	y[1:4] = g(x[1:4])
	for (i in 4:(n-1))
		y[i+1] <- it.adams(h, c(x[i-3], x[i-2], x[i-1], x[i], x[i+1]), c(y[i-3], y[i-2], y[i-1], y[i]))
	
	#par(mfcol = c(1, 2))
	plot(x, y, type="l", col="blue")
	title("Adams")
	#lines(x, g(x), type="l", col="blue")
	for (i in seq(from=1, to=length(x), by=floor(length(x)/10)))
		mprint(x[i], y[i])
}

it.adams <- function(h, x, y){
	res <- y[4] + h * (55/24 * f(x[4], y[4]) - 59/24 * f(x[3], y[3]) + 37/24 * f(x[2], y[2]) - 9/24 * f(x[1], y[1]))
	
	return ( y[4] + h * (9/24 * f(x[5], res) + 19/24 * f(x[4], y[4]) - 5/24 * f(x[3], y[3]) + 1/24 * f(x[2], y[2])) )
}

mprint <- function(x, y){
	color <- "red"
	if (abs(g(x) - y) < .001) color<-"blue"
	if (abs(g(x) - y) < .00001) color<-"green"
	points(x, y, pch=17, col=color)
	#print(abs(g(x) - y))
}

runge.kutta.adams.diff <- function(a, b, n.c){
	par(mfcol = c(1, 2))
	runge.kutta(a, b, n.c, FALSE)
	adams(a, b, n.c)
}