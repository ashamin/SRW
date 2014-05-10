# TODO: Add comment
# 
# Author: Антон
###############################################################################

data <- data.frame(x=x, y=y)
rm(x, y)
attach(data)

fx <- function(){
	detach(data)
	fix(data)
	#always sort. can't overwhelm it
	#data <- data.frame(x=as.real(levels(factor(data$x))), y=as.real(levels(factor(data$y))))
	dup <- factor(data$x)
	dataY <- split(data$y, dup)
	y <- numeric(length(dataY))
	x <- numeric(length(dataY))
	for (i in 1:length(dataY)){
		x[i] <- as.real(levels(dup)[i])
		y[i] <- dataY[[i]][1]
	}
	data <- data.frame(x=x, y=y)
	attach(data)
}