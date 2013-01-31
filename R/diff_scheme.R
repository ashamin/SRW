# TODO: Solves equation 
#
#			u'' + A(x)*u' + B(x)*u = C(x)
#
#			with border conditions:
#
#				F1*u(a) + D1*u'(a) = E1
#				F2*u(b) + D2*u'(b) = E2
#
#			x from interval [a, b]
# 
# Author: ?????
###############################################################################

#derivative equation coefficients A(x), B(x), C(x)
A <- function(x){
	-2
}

B <- function(x){
	0
}

C <- function(x){
	exp(x)*(x^2 + x - 3)
}

#coefficients of border conditions
a <- 0
b <- 1
F1 <- 0
D1 <- 1
E1 <- 2
F2 <- -1
D2 <- 1
E2 <- exp(1) * (exp(1) - 3)

#solves matrix with 3 diagonals with sweep method
#@parameters a - lower diag, b - diag, c - upper diag, d - vector
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

#forming full scheme without border conditions
form.scheme <- function(d1, d2, d3, x, h){
	n <- length(x)
	
	#lower diag
	d1[1:(n-2)] <- (1/h  -  A(x[2:(n-1)])/2) / h
	#main diag
  d2[2:(n-1)] <- (B(x[2:(n-1)])  -  2/h^2)
	#upper diag
  d3[2:(n-1)] <- (1/h  +  A(x[2:(n-1)])/2) / h
	
	return (list(d1=d1, d2=d2, d3=d3))
}

#main method of forming differential scheme
diff_method <- function(n.c, precision){
	
	if (precision != 1 && precision != 2){
		print(paste("no formula for ", as.character(precision), " precision"))
		return()
	}
	
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	d1 <- c(rep(0, times=n.c))
	d2 <- c(rep(0, times=n.c+1))
	d3 <- c(rep(0, times=n.c))
	f <- C(x)
	
	if (precision == 1){
		f[1] <- E1
		f[n.c+1] <- E2
		d2[1] <- F1 - D1/h
		d3[1] <- D1/h
		d1[n.c] <- -D2/h
		d2[n.c+1] <- F2 + D2/h
	}
	else if (precision == 2){
		f[1] <- E1 + h*D1*C(a)/2
		f[n.c+1] <- E2 - h*D2*C(b)/2
		d2[1] <- F1 - D1/h + B(a)*D1*h/2 - A(a)*D1/2
		d3[1] <- D1/h + A(a)*D1/2
		d1[n.c] <- -D2/h + D2*A(b)/2
		d2[n.c+1] <- F2 + D2/h - B(b)*D2*h/2 -A(b)*D2/2 
	}
	
	diags <- form.scheme(d1, d2, d3, x, h)
	
	
	
	
	
	#print(diags)
	#print(f)
	#print(sweep(diags$d3, diags$d2, diags$d1, f))
	

	M <- matrix(rep(0, times=length(x)^2), length(x), length(x))
	
	dd1 <- t(matrix(rep(1:length(M[1,]), each=2), 2, length(M[1,])))
	dd2 <- matrix(c(1:(length(x)-1), 2:length(x)), length(x)-1, 2)
	dd3 <- matrix(c(2:length(x), 1:(length(x)-1)), length(x)-1, 2)
	
	M[dd1] <- diags$d2
	M[dd2] <- diags$d3
	M[dd3] <- diags$d1
  
  source("R/TDMA.R")
  
	#return (solve(M, f))
  return (TDMA(diags$d1, diags$d2, diags$d3, f))
}

test <- function(){
	n.c <- 25
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	n.c <- 26
	res <- exp(x)*(exp(x) - x^2 - x + 1)
	p.ans <- C(x[2:n.c])
	a1.ans <- c(rep(0, times=n.c-1))
	#for (i in 2:2)
	#a1.ans[i-1] <- (res[i+1] - 2*res[i] + res[i-1])/h^2 + A(x[i])*(res[i+1] - res[i-1])/(2*h) + B(x[i])*res[i]
	a1.ans <- (res[2:(n.c-1)+1] - 2*res[2:(n.c-1)] + res[2:(n.c-1)-1]) / h^2 + A( x[2:(n.c-1)])*(res[2:(n.c-1)+1] - res[2:(n.c-1)-1]) / (2*h) +  B(x[2:(n.c-1)]) * res[2:(n.c-1)]
	a2.ans <- res[2:(n.c-1)-1] * (1/h - A( x[2:(n.c-1)])/2)/h   +   res[2:(n.c-1)]*(B(x[2:(n.c-1)]) - 2/h^2)   +   res[2:(n.c-1)+1]*(1/h + A( x[2:(n.c-1)])/2)/h
	dgs <- diff_method(100, 1)
	a3.ans <- res[2:(n.c-1)-1] * (dgs$d1[1:(n.c-2)])  +  res[2:(n.c-1)] * dgs$d2[2:(n.c-1)]  +  res[2:(n.c-1)+1] * dgs$d3[2:(n.c-1)]
	print(list(res=res, p.ans=p.ans, a1.ans=a1.ans, a2.ans=a2.ans, a3.ans=a3.ans))

}

#build two diff schemes with different precision
# @param n.c number of intervals in differential scheme
diff <- function(n.c){
	h <- (b-a)/n.c
	x <- seq(from=a, to=b, by=h)
	p1 <- diff_method(n.c, 1)
	p2 <- diff_method(n.c, 2)
	xx <- seq(from=a, to=b, by=.025)
	par(mfcol = c(1, 2))
	plot(c(0,1), c(2, 5), pch=17)
	lines(xx, exp(xx)*(exp(xx) - xx^2 - xx + 1))
	for (i in 1:length(x))
		mprint(p1[i], x[i])
	plot(c(0,1), c(2, 5), pch=17)
	lines(xx, exp(xx)*(exp(xx) - xx^2 - xx + 1))
	for (i in 1:length(x))
		mprint(p2[i], x[i])
}



tst <- function(){
	M <- matrix(rep(0, 25), 5, 5)
	for (i in 1:5)
		M[i, i] <- as.double(19)
	M[1,2] <- 2
	M[2,3] <- 3
	M[3,4] <- 8
	M[4,5] <- 1
	
	M[2, 1] <- 0.6
	M[3, 2] <- 1
	M[4, 3] <- 1.1
	M[5,4] <- 7
	
	print(M)
	
	b <- c(1, 4, 5, 2, 1)
	
	print(solve(M, b))
	
	print(sweep(c(0.6, 1, 1.1, 7), c(rep(19, times=5)) ,c(2, 3, 8, 1), b))
}


mprint <- function(res1, x){
	color <- "red"
	if (abs(res1 - exp(x)*(exp(x) - x^2 - x + 1)) < .001) color<-"blue"
	if (abs(res1 - exp(x)*(exp(x) - x^2 - x + 1)) < .00001) color<-"green"
	points(x, res1, pch=17, col=color)
}
