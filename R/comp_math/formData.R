# TODO: Add comment
# 
# Author: Антон
###############################################################################

f <- function(x0){
	abs(x0)
}

x <- seq(from=-4, to=4, by=.5)#c(seq(from=-1, to=0, by=.1), .5, 1)#c(-1, -.25, 0.25, 1, -.1, .1)
y <- f(x)

spline <- data.frame(a=real(length(x)), b=real(length(x)), c=real(length(x)), d=real(length(x)))