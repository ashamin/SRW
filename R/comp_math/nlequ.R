# TODO: Numeric Analisys. Work 11. Solving nonlinear equations
# 
# Author: Антон
###############################################################################

#epsilon <- .00000001
#step <- .0000001

epsilon <- .0000001
step <- .000001

f <- function(x0){
	#-x0^2 + x0 - 2
	#x0^3 + 2*x0^2 - 5*x0 - 6
	#x0^3
	#-(x0^(-2))+2
	#(x0-2)^(-2)-5
	#x0^2
	#2*(x0+3)^2-8
	#-x0^6 + (x0+1)^3 + (x0-2)^2 +5
	#-x0^2 + 2
	#-sin(x0)^2+.3
	#sin((x0+1)^2)/3
	sin(x0^7)*x0^3+cos(x0^2)
	#sin(x0^2)
	#tan(x0+3)/2
	#atan(x0+2)^2-1
}

rdev1 <- function(f, x0){
	rh <- .000001
	(f(x0+rh/2) - f(x0-rh/2))/ rh
}

rdev2 <- function(f, x0){
	rh <- .000001
	(f(x0+rh) - 2*f(x0) + f(x0-rh)) / (rh^2)
}

sign <- function(x0){
	if (x0 > 0) return(1)
	else return(0)
}

draw <- function(a, b){
	x <- seq(from=a, to=b, by=.025)
	plot(c(a, b), c(a, b), xlab="X axis", ylab="Y axis")
	lines(x, rep(0, times=length(x)))
	lines(x, f(x))
}

#splits interval [a, b] to subintervals and search roots on every subinterval 
# add if we found root on a border i++ or i+2
ispg <- function(a, b){
	draw(a, b)
	rootnumber <<- 1
	x <- seq(from=a, to=b, by=step)
	#for (i in 1:(length(x)-1)){
	i <- 1
	while (i < length(x)){
		if (f(x[i]) == 0){
			mprint(x[i])
			i <- i+1
		} 
		else if (f(x[i+1]) == 0){
			mprint(x[i+1])
			i <- i+2
		} 
		else if (sign(f(x[i])) != sign(f(x[i+1]))) mprint(comb(x[i], x[i+1]))
		i <- i+1
	}
	print_roots()
}

print_roots <- function(){
	roots <- data.frame(root=rts, residual=res)
	rts <<- double(0)
	res <<- double(0)
	print(roots)
}

rts <- double(0)
res <- double(0)

mprint <- function(x){
	color <- "red"
	if (abs(f(x)) < .001) color<-"blue"
	if (abs(f(x)) < .00001) color<-"green"
	points(x, f(x), pch=17, col=color)
	
	rts <<- c(rts, x)
	res <<- c(res, f(x))
}

#add options with different combinations of signs of f'(x) and f''(x)
#combined chord and tangent(Newton) methods
comb <- function(a, b){
	xN <- b
	xC <- a
	if (rdev2(f, xN)*rdev1(f, xN) > 0){
		while((xN - xC) >= epsilon){
			xC <- xC - f(xC)*(xN-xC)/(f(xN) - f(xC))
			xN <- xN - f(xN)/rdev1(f, xN) 
		}
	}
	else{
		while((xN - xC) >= epsilon){
			xN <- xN - f(xN)*(xC-xN)/(f(xC) - f(xN))
			xC <- xC - f(xC)/rdev1(f, xC) 
		}
	}
	xN
}

