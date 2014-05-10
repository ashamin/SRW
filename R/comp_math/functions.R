# TODO: Numeric Analysis. Works 6,7. Lagrange polynomial & spline interpolation
# 		Use after formData.R
#		Can be used in cooperate with dataFrame.R. Functions from that file provides
#			edit coordinates of interpolation nodes.
#
# Author: ?????
###############################################################################

aa <- function(x0){
	list(interpolation=int(x0), spline_value=spl(x0), real_value=f(x0))
}

int <- function(x0){
	if (!is.vector(x0)){
		L <- real(length(x))
		for(k in 1:length(x))
			L[k] <- y[k] * prod(-x[-k] + x0) / prod(-x[-k]+x[k])
		sum(L)
	}
	else{
		L <- real(length(x))
		v <- real(length(x0))
		for (i in 1:length(x0)){
			for(k in 1:length(x))
				L[k] <- y[k] * prod(-x[-k] + x0[i]) / prod(-x[-k]+x[k])
			v[i] <- sum(L)
		}
		v
	}
}

#solves matrix with 3 diagonals with sweep method
#@parameters a - upper diag, b - diag, c - lower diag, d - vector
sweep <- function(a, b, c, d){
	n <- length(d)
	x <- real(n)
	alpha <- c(rep(0, times=n+1)) 
	betta <- c(rep(0, times=n+1))
	a <- c(0, a)
	c <- c(c, 0)
	alpha[2:(n+1)] <- c[1:n] / (b[1:n] - a[1:n]*alpha[1:n])
	betta[2:(n+1)] <- (a[1:n]*betta[1:n]-d[1:n]) / (b[1:n]-a[1:n]*alpha[1:n])
	
	x[n] <- betta[n+1]
			
	x[(n-1):1] <- alpha[n:2]*x[n:2] - betta[n:2]
	x
}

#compute spline functions
cspl <- function(){
	attach(spline)
	a <- c(y)
	h <- c(x[2:length(x)]-x[1:(length(x)-1)])
	
	c <- real(length(x)-2)
	M <- matrix(rep(0, times=length(c)^2), length(c), length(c))
	
	d1 <- t(matrix(rep(1:length(M[1,]), each=2), 2, length(M[1,])))
	d2 <- t(matrix(c(1:(length(c)-1), 2:length(c)), 2, length(c)-1))
	d3 <- matrix(c(2:length(c), 1:(length(c)-1)), length(c)-1, 2)
	
	M[d2] <- h[2:length(c)]
	M[d1] <- 2*(h[1:length(c)]+h[2:length(h)])
	M[d3] <- h[2:length(c)]
	
	r <- 6*(((y[3:length(x)]-y[2:(length(x)-1)]) / h[2:length(h)]) - ((y[2:(length(x)-1)]-y[1:(length(x)-2)]) / h[1:(length(h)-1)]))
	#c <- c(0, solve(M,r), 0)
	c <- c(0, sweep(M[d2], M[d1], M[d3], r), 0)
	
	d <- (c[2:length(c)] - c[1:(length(c)-1)]) / h[1:length(h)]
	b <- h[1:length(h)]*c[2:length(c)] / 2 - h[1:length(h)]^2 * d[1:length(d)] / 6 + (y[2:length(y)] - y[1:(length(y)-1)]) / h[1:length(h)]
	
	spline <- data.frame(a=a[2:length(a)], b=b, c=c[2:length(c)], d=d)
	detach()
	spline
}

spline <- cspl()

spl <- function(x0){
	for (i in 2:length(spline$a))
		if (x[i-1] <= x0 & x[i] >= x0){
			i <- i-1
			break
		}
	 
	 (spline$a[i] + spline$b[i]*(x0 - x[i+1]) + spline$c[i]*((x0 - x[i+1])^2)/2 + spline$d[i]*((x0 - x[i+1])^3)/6)
}

draw <- function(){
	xl <- seq(from=min(x), to=max(x), by=.025)
	plot(c(x, min(y), max(y)), c(y, min(x), max(x)))
	lines(xl, f(xl), col="black")
	lines(xl, int(xl), col="red")
	spline <<- cspl()
	for (i in 2:length(xl))
		lines(c(xl[i-1], xl[i]), c(spl(xl[i-1]), spl(xl[i])), col='green')
}
