# TODO: Numeric Analysis. Work 13. Nonlinear equation systems
#			!Convergence of method strongly depends from initial приближения
# 
# Author: Антон
###############################################################################

f1 <- function(x){
	x[1]^2 + x[2]^2 + x[3]^2 - 1
}

f2 <- function(x){
	2*x[1]^2 + x[2]^2 - 4*x[3]
}

f3 <- function(x){
	3*x[1]^2 - 4*x[2] - x[3]^2
}

F <- function(x){
	c(f1(x), f2(x), f3(x))
}

j11 <- function(x){
	2*x[1]
}
j12 <- function(x){
	2*x[2]
}
j13 <- function(x){
	2*x[3]
}

j21 <- function(x){
	4*x[1]
}
j22 <- function(x){
	2*x[2]
}
j23 <- function(x){
	-4
}

j31 <- function(x){
	6*x[1]^2
}
j32 <- function(x){
	-4
}
j33 <- function(x){
	-2*x[3]
}

J <- function(x){
	t(matrix(c(j11(x), j12(x), j13(x), j21(x), j22(x), j23(x), j31(x), j32(x), j33(x)), 3, 3))
}

NewtonSolve <- function(x0=rnorm(3), epsilon=rep(.0001, times=3)){
	rts <- it.process(x0, epsilon)
	res <- F(rts)
	message("...Roots & Residuals")
	print(data.frame(x0=x0, I="|", root=rts, i="|", residual=res))
}

it.process <- function(x0, epsilon, maxit=1000){
	x <- x0
	it <- 0
	repeat{
		it <- it+1
		if (it > maxit){
			print(data.frame(x0=x0))
			stop("iteration process obviously won't converge")
		}
		xp <- solve(J(x), -F(x))
		x <- x + xp
		if (prod(abs(xp) < epsilon) == 1) break;
	} 
	return(x)
}
