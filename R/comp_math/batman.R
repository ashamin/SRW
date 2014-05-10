# TODO: Numerical Analysis. Work 10. Presetnation of geometrical Monte-Carlo method
# 
# Author: Антон
###############################################################################


f1 <- function(x){
	2*sqrt((-abs(abs(x)-1))*abs(3-abs(x))/((abs(x)-1)*(3-abs(x))))*(1+abs(abs(x)-3)/(abs(x)-3))*sqrt(1-(x/7)^2)+(5+.97*(abs(x-.5)+abs(x+.5))-3*(abs(x-.75)+abs(x+.75)))*(1+abs(1-abs(x))/(1-abs(x)))
}
f2 <- function(x){
	(2.71052+1.5-.5*abs(x)-1.35526*sqrt(4-(abs(x)-1)^2))*sqrt(abs(abs(x)-1)/(abs(x)-1))+.9
}
f3 <- function(x){
	(-3)*sqrt(1-(x/7)^2)*sqrt(abs(abs(x)-4)/(abs(x)-4))
}
f4 <- function(x){
	abs(x/2)-.0913722*x^2-3+sqrt(1-(abs(abs(x)-2)-1)^2)
}

batman <- function(){
	plot(c(-8, 8), c(-8,8))
	x <- seq(from=-7, to=-4.025, by=.025)
	y <- f3(x)
	lines(x, y)
	x <- seq(from=4.025, to=7, by=.025)
	y <- f3(x)
	lines(x, y)
	x <- seq(from=-4, to=4, by=.025)
	y <- f4(x)
	lines(x, y)
	x <- seq(from=-3, to=-1.025, by=.025)
	y <- f2(x)
	lines(x, y)
	x <- seq(from=1.025, to=3, by=.025)
	y <- f2(x)
	lines(x, y)
	x <- seq(from=-7, to=-3.025, by=.025)
	y <- f1(x)
	lines(x, y)
	x <- seq(from=-.99, to=.99, by=.01)
	y <- f1(x)
	lines(x, y)
	x <- seq(from=3.025, to=7, by=.025)
	y <- f1(x)
	lines(x, y)
	S()
}

putdot <- function(x, y, color){
	for (i in 1:length(x))
		lines(c(x[i], x[i]+.001), c(y[i], y[i]+.001), col=color)
}

#premium for split with cut demostration
S1 <- function(){
	K <- 0
	N <- 30000
	x <- runif(N, -7, 7)
	y <- runif(N, -4, 4)
	sx <- split(x, cut(x, c(-7, -4, -3, -1, 1, 3, 4, 7)))
	sy <- split(y, cut(x, c(-7, -4, -3, -1, 1, 3, 4, 7)))
	for(i in 1:length(sx$'(-7,-4]')){
		if (sy$'(-7,-4]'[i] >= f3(sx$'(-7,-4]'[i]) && sy$'(-7,-4]'[i] <= f1(sx$'(-7,-4]'[i])){
			putdot(sx$'(-7,-4]'[i], sy$'(-7,-4]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(-7,-4]'[i], sy$'(-7,-4]'[i], "black")
	}
	for(i in 1:length(sx$'(-4,-3]')){
		if (sy$'(-4,-3]'[i] >= f4(sx$'(-4,-3]'[i]) && sy$'(-4,-3]'[i] <= f1(sx$'(-4,-3]'[i])){
			putdot(sx$'(-4,-3]'[i], sy$'(-4,-3]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(-4,-3]'[i], sy$'(-4,-3]'[i], "black")
	}
	for(i in 1:length(sx$'(-3,-1]')){
		if (sy$'(-3,-1]'[i] >= f4(sx$'(-3,-1]'[i]) && sy$'(-3,-1]'[i] <= f2(sx$'(-3,-1]'[i])){
			putdot(sx$'(-3,-1]'[i], sy$'(-3,-1]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(-3,-1]'[i], sy$'(-3,-1]'[i], "black")
	}
	for(i in 1:length(sx$'(-1,1]')){
		if (sy$'(-1,1]'[i] >= f4(sx$'(-1,1]'[i]) && sy$'(-1,1]'[i] <= f1(sx$'(-1,1]'[i])){
			putdot(sx$'(-1,1]'[i], sy$'(-1,1]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(-1,1]'[i], sy$'(-1,1]'[i], "black")
	}		
	for(i in 1:length(sx$'(1,3]')){
		if (sy$'(1,3]'[i] >= f4(sx$'(1,3]'[i]) && sy$'(1,3]'[i] <= f2(sx$'(1,3]'[i])){
			putdot(sx$'(1,3]'[i], sy$'(1,3]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(1,3]'[i], sy$'(1,3]'[i], "black")
	}
	for(i in 1:length(sx$'(3,4]')){
		if (sy$'(3,4]'[i] >= f4(sx$'(3,4]'[i]) && sy$'(3,4]'[i] <= f1(sx$'(3,4]'[i])){
			putdot(sx$'(3,4]'[i], sy$'(3,4]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(3,4]'[i], sy$'(3,4]'[i], "black")
	}
	for(i in 1:length(sx$'(4,7]')){
		if (sy$'(4,7]'[i] >= f3(sx$'(4,7]'[i]) && sy$'(4,7]'[i] <= f1(sx$'(4,7]'[i])){
			putdot(sx$'(4,7]'[i], sy$'(4,7]'[i], "yellow")
			K <- K+1
		}
		else putdot(sx$'(4,7]'[i], sy$'(4,7]'[i], "black")
	}
	print(14*8*K/N)
}

S <- function(){
	N <- 30000
	K <- 0
	x <- runif(N, -7, 7)
	y <- runif(N, -4, 4)
	
	for(i in 1:length(x)){
		
		if (x[i] > -7 && x[i]< -4){
			if (y[i] >= f3(x[i]) && y[i] <= f1(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		}
		
		if (x[i] > -4 && x[i]< -3) 
			if (y[i] >= f4(x[i]) && y[i] <= f1(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		
		if (x[i] > -3 && x[i] < -1)
			if (y[i] >= f4(x[i]) && y[i] <= f2(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		
		if (x[i] > -1 && x[i] < 1)
			if (y[i] >= f4(x[i]) && y[i] <= f1(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		
		
		if (x[i] > 1 && x[i] < 3)
			if (y[i] >= f4(x[i]) && y[i] <= f2(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		
		
		if (x[i] > 3 && x[i] < 4) 
			if (y[i] >= f4(x[i]) && y[i] <= f1(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
		
		if (x[i] > 4 && x[i]<7)
			if (y[i] >= f3(x[i]) && y[i] <= f1(x[i])){
				putdot(x[i], y[i], "yellow")
				K <- K+1
			}
			else putdot(x[i], y[i], "black")
	}
	print(14*8*K/N)
}