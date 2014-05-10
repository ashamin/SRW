# TODO: Numeric analysis. Work 12. Muller && Parab methods
# 
# Author: ¿ÌÚÓÌ
###############################################################################
#source("Lagrange polynomial/parab.R")
epsilon <- .0001

#initfunc <- c(1, 2, -5, -6)
#initfunc <- c(1, -2)
#initfunc <- c(2, 5, 4, 1, -7)
initfunc <- c(.5462, .987654, .1235, -.234547, -.6762, .657885, -.324235, -.212434, -.776661, -.2457684)
func <- initfunc

f <- function(x0){
	sum(func[1:length(func)]*x0^((length(func):1)-1))+0i
}

parab <- function(){
	func <<- initfunc
	rts <- complex(length(initfunc)-1)
	its <- integer(length(initfunc)-1)
	for (i in 1:(length(initfunc)-2)){
		xn <- 0#-2 + 0i
		xn1 <- 1#-3 + 0i
		xn2 <- -1#5 + 0i
		it <- 0
		while(diff(xn, xn1) >= epsilon && it < 1000 && abs(f(xn)) > epsilon){
			it <- it +1
			#compute coefs & roots
			cfs <- coefs(c(xn, xn1, xn2))
			a <- cfs$a
			b <- cfs$b
			cc <- cfs$cc
			#attach(cfs)
			xn2 <- xn1 
			xn1 <- xn
			############ «ƒ≈—‹ ¬—≈√ƒ¿ œ–»¡¿¬ÀﬂÀ» ÃŒƒ”À‹!  Õ”∆ÕŒ œ–»¡¿¬Àﬂ“‹ ¿¡—ŒÀﬁ“ÕŒ≈ «Õ¿◊≈Õ»≈
			if ( mn((-b+sqrt(b^2 - 4*a*cc))/(2*a), (-b-sqrt(b^2 - 4*a*cc))/(2*a)) == abs((-b+sqrt(b^2 - 4*a*cc))/(2*a)) )
				xn <- xn1 + (-b+sqrt(b^2 - 4*a*cc))/(2*a)
			else xn <- xn1 + (-b-sqrt(b^2 - 4*a*cc))/(2*a)
			#detach(cfs)
		}
		if (it == 1000) print("D'Oh")
		rts[i] <- xn
		its[i] <- it
		rm_root(xn)
	}
	rts[length(rts)] <- -func[2]/func[1]
	its[length(its)] <- 1
	func <<- initfunc
	res <- complex(length(rts))
	for (i in 1:length(res))
		res[i] <- f(rts[i])
	data.frame(root=rts, I = rep("|", times=length(rts)), iterations=its, i = rep("|", times=length(rts)), residual=res)
}

Muller <- function(){
	func <<- initfunc
	rts <- complex(length(initfunc)-1)
	its <- integer(length(initfunc)-1)
	for (i in 1:(length(initfunc)-2)){
		xn <- 0#-2 + 0i
		xn1 <- 1#-3 + 0i
		xn2 <- -1#5 + 0i
		it <- 0
		while(diff(xn, xn1) >= epsilon && it < 1000 && abs(f(xn)) > epsilon){
			it <- it +1
			#compute coefs & roots
			cfs <- coefs(c(xn, xn1, xn2))
			a <- cfs$a
			b <- cfs$b
			cc <- cfs$cc
			#attach(cfs)
			xn2 <- xn1 
			xn1 <- xn
			############ «ƒ≈—‹ ¬—≈√ƒ¿ œ–»¡¿¬ÀﬂÀ» ÃŒƒ”À‹!  Õ”∆ÕŒ œ–»¡¿¬Àﬂ“‹ ¿¡—ŒÀﬁ“ÕŒ≈ «Õ¿◊≈Õ»≈
			if ( mx(b + sqrt(b^2 - 4*a*cc), b - sqrt(b^2 - 4*a*cc)) == abs(b + sqrt(b^2 - 4*a*cc)) )
				xn <- xn1 - 2*cc / ( b + sqrt(b^2 - 4*a*cc) )
			else xn <- xn1 - 2*cc / ( b - sqrt(b^2 - 4*a*cc) )
			#detach(cfs)
		}
		if (it == 1000) print("D'Oh")
		rts[i] <- xn
		its[i] <- it
		rm_root(xn)
	}
	rts[length(rts)] <- -func[2]/func[1]
	its[length(its)] <- 1
	func <<- initfunc
	res <- complex(length(rts))
	for (i in 1:length(res))
		res[i] <- f(rts[i])
	data.frame(root=rts, I = rep("|", times=length(rts)), iterations=its, i = rep("|", times=length(rts)), residual=res)
}

pM <- function(){
	digs <- 6
	p <- parab()
	M <- Muller()
	if (diffs(p, M) != 0)
		print(M, digits=digs)
	print(p, digits=digs)
}

diffs <- function(p, M){
	print("Differences between parab method & Muller method")
	pr <- 0
	for (i in 1:length(p$root))
		if (diff(p$root[i], M$root[i]) > .000001){
			print(paste(as.character(p$root[i]), "|", as.character(M$root[i])))
			pr <- 1
		} 
	if (pr == 0) print("None")
	pr
}

mx <- function(x1, x2){
	max(sqrt(Re(x1)^2 + Im(x1)^2), sqrt(Re(x2)^2 + Im(x2)^2))
}

mn <- function(x1, x2){
	min(sqrt(Re(x1)^2 + Im(x1)^2), sqrt(Re(x2)^2 + Im(x2)^2))
}

diff <- function(x1, x2){
	#abs(sqrt(Re(x1)^2 + Im(x1)^2) - sqrt(Re(x2)^2 + Im(x2)^2))
	abs(Re(x1) - Re(x2)) + abs(Im(x1) - Im(x2))
}

#computes coefficient of Newton polynom. x=(xn, xn-1, xn-2)
#returns coefficients list (a=square, b=linear, c=free)
coefs <- function(x){
	list(a=div_diff(x), b=sum(div_diff(x[-3]), div_diff(x[-2]), -div_diff(x[-1])), cc=f(x[1]))
	#list(a=div_diff(x), b=div_diff(x)*(x[1] - x[2])+div_diff(x[-3]), c=f(x[1]))
}

#removes root from function
rm_root <- function(root){
	res <- complex(length(func)-1)
	for (i in 1:(length(func)-1)){
		res[i] <- func[i]
		func[i+1] <- func[i+1] + root*func[i]
		func[i] <- 0
	}
	func <<- res
}

#computes divided differece of vector x
div_diff <- function(x){
	sum <- 0 + 0i
	for (k in 1:length(x))
		sum <- sum + f(x[k])/prod(x[k]-x[-k])
	return (sum)
}

