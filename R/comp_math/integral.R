# TODO: Numeric Analysis. Work 9. 
# 		R can't pass recursion big depth
# Author: Антон
###############################################################################

#options(expr = 50000)

f <- function(x0){
	#sin(x0)
	#x0^9+x0^3/8
	x0^2
}

hs <- 0.1
epsilon = .0001
steps <- 1

sq <- function(a, b, h=hs){
	x <- seq(from=a, to=b, by=h)
	sum(abs(f(x[(2:length(x))]-h/2)))*h
}

tr <- function(a, b, h=hs){
	x <- seq(from=a, to=b, by=h)
	(sum(abs(f(x[c(1, length(x))])/2)) + sum(abs(f(x[-c(1, length(x))])))) * h
}

Simpson <- function(a, b, h=hs){
	x <- seq(from=a, to=b, by=h)
	n <- length(x)
	sum(abs(f(x[1:(n-1)])) + abs(f(x[2:n])) + 4*abs(f(x[2:n]-h/2)))*h/6
}

#auto computing step of integration
#m - в формуле погрешности - степень h в порядке точности формулы M
auto <- function(funct, a, b, m, h=(b-a)){
	x <- seq(from=a, to=b, by=h)
	n <- length(x)
	res <- real(n-1)
	diff <- real(n-1)
	for(i in 1:(n-1))
		res[i] <- funct(x[i], x[i+1], h)
	
	for(i in 1:(n-1))
		#diff[i] <- (abs(funct(x[i], x[i+1], h/2) - res[i])*(b-a)) / ((2^m-1)*h)
		diff[i] <- abs(funct(x[i], x[i+1], h/2) - res[i])
	
	for(i in 1:(n-1))
		if (diff[i] > epsilon){
			steps <<- steps+1
			if (h > .0001) res[i] <-auto(funct, x[i], x[i+1], m, h/2)
		}
	sum(abs(res))
}

comp <- function(a, b){
	res <- double(6)
	st <- integer(6)
	res[1]  <- sq(a, b)
	st[1] <- (b-a)/hs
	steps <<- 1
	res[2] <- auto(sq, a, b, 2)
	st[2] <- steps
	res[3]  <- tr(a, b)
	st[3] <- (b-a)/hs
	steps <<- 1
	res[4] <- auto(tr, a, b, 2)
	st[4] <- steps
	res[5]  <- Simpson(a, b)
	st[5] <- (b-a)/hs
	steps <<- 1
	res[6] <- auto(Simpson, a, b, 4)
	st[6] <- steps	
	names(res)[1] <- "sq"
	names(res)[2] <- "auto sq"
	names(res)[3] <- "tr"
	names(res)[4] <- "auto tr"
	names(res)[5] <- "Simpson"
	names(res)[6] <- "auto Simpson"
	names(st)[1] <- "sq"
	names(st)[2] <- "auto sq"
	names(st)[3] <- "tr"
	names(st)[4] <- "auto tr"
	names(st)[5] <- "Simpson"
	names(st)[6] <- "auto Simpson"
	print(t(res))
	print(t(st))
}
