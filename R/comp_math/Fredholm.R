# TODO: Add comment
# 
# Author: Антон
###############################################################################

k <- function(t1, t2){
	t1 * t2
}

f <- function(t){
	2 * t / 3
}

uSolve <- function(t){
	t	
}

Simpson.coef <- function(h, n){
	coefs <- rep(0, times=n)
	coefs[-(2:(n-1))] = h/3
	coefs[seq(from=2, to=n-1, by=2)] = 4*h/3
	coefs[seq(from=3, to=n-1, by=2)] = 2*h/3
	return (coefs)
}

tr.coef <- function(h, n){
	coefs <- rep(0, times=n)
	coefs[-(2:(n-1))] = h/2
	coefs[2:(n-1)] = h
	return (coefs)
}

sq.coef <- function(h, n){
	return (rep(h, times=n))
}

# Solves a Fredhold integral equation of II type
#
# Make and solve linear system of equaions with result corresponds to values of 
# u(x) function in node points. Node points computes automatically using [a, b]
# step value h. 
#
# @param a left interval border
# @param b right interval border
# @param h step value
# @param lambda labda value
# @param K K function
# @param f f function
# @param coef.method method of computing quadratic equation for calculation integrals
#
#
# @return values of u(x) function in node points
#
solve.Fredholm <- function(a, b, h, lambda, K=k, fr=f, coef.method=tr.coef){
	t <- seq(from=a, to=b, by=h)
	n <- length(t)
	coefs <- coef.method(h, n)
	A <- matrix(coefs, n, n, byrow=TRUE)
	A <- A * lambda;
	
	for (i in 1:n)
		A[i, 1:n] <- A[i, 1:n] * K(t[i], t[1:n]) * (-1)
	
	#A[1:n, 1:n] <- A[1:n, 1:n] * K(t[1:n], t[1:n]) * (-1)
	#K(t[1:n], t[1:n]) ------ считает не то! потому что функция, а не массив
	
	for (i in 1:n)
		A[i, i] <- 1 + A[i, i]
	
	F <- fr(t)
	result = solve(A, F);
	draw(t, result)
	list(A=A, b=F, result=solve(A, F), residuals.norm=comp.residuals(t, result), coefs=coefs)
}

draw <- function(t, res){
	spl <- splinefun(t, res)
	#plot(t, res, xlim=c(0, 3), ylim=c(0,3), pch=17, col="red")
	#x <- seq(from=min(t), to=max(t), by=.025)
	#lines(x, spl(x))
	plot(t, res, xlim=c(0, 3), ylim=c(0,3))
	curve(spl)
}

#
# Compute norm of solve error.
#
comp.residuals <- function(t, result){
	return (max(abs(result - uSolve(t))))
}

Fredholm <- function(a, b, h, lambda, K=k, fr=f, coef.method=tr.coef){
	residuals <- data.frame(SQ = solve.Fredholm(a, b, h, lambda, K, fr, sq.coef)$residuals, i="|",
						TR = solve.Fredholm(a, b, h, lambda, K, fr, tr.coef)$residuals, I="|",
							SIMPSON = solve.Fredholm(a, b, h, lambda, K, fr, Simpson.coef)$residuals)
	print (residuals)
}
