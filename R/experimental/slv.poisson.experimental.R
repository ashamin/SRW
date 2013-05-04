


print.temperaturegrams.par.method <- function(x, y, A, b, I, J, n, h){
  ra <- matrix(rep(0, (length(x)-2)*(length(y)-2)), length(x)-2, length(y)-2, byrow=TRUE)
  for (i in 2:(length(x)-1))
    for (j in 2:(length(y)-1))
      ra[i-1, j-1] = right.answer(x[i], y[j])
  par(mfcol=c(1, 2))
  
#   P=SSOR.par(1, n, I, h)
#   print(P$Mx)
#   print(P$My)
  
  s <- mincorr.par(A, b, P=SSOR.par(.5, n, h), maxit=10000)
  
  #   s <- mincorr(A, b, P=SSOR.precond.approach1(A, 1), maxit=1000))
  #   s <- minres(A, b, P=SSOR.precond.approach1(A, 1), maxit=10000)
  #   s <- minres(A, b, P=SSOR.precond.approach1(A, 1), maxit=1000, tau.compute=FALSE)
  
  sl <- matrix(s$root, (I-2), (J-2), byrow=FALSE)
  image(x[2:(I-1)], y[2:(J-1)], sl)
  title("solve")
  image(x[2:(I-1)], y[2:(J-1)], ra)
  title("right answer")
  
  print(data.frame(residual.norm=norm(ra - sl)/norm(sl), iteration.number=unique(s$it)))  
}