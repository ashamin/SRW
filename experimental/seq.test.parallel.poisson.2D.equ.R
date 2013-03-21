
# works with version of solving poisson 2D equation in R/poisson.equ.sq.2D.R implementation

source("R/poisson.equ.sq.2D.R")

# test answer function
# answ <- function(x, y){
#   -x^3/6
# }

slv.seq.par <- function(a=a.test, b=b.test, g1=g1.test, g2=g2.test,
                        g3=g3.test, g4=g4.test, f=f.test,
                        I=10, J=10, h1=a/(I-1), h2=b/(J-1)){
                        
  x <- seq(from=0, to=a, by=h1)
  y <- seq(from=0, to=b, by=h2)
  
  n <- I-2
  m <- J-2
  
  res1 <- matrix(rep(0, n*m), n, m) 
  res2 <- matrix(rep(0, n*m), n, m) 
  
#   A <- form.A.matrix.1(n, h1, I)
#   b <- form.b.1(n, I, x, y[2], h1, function(x, y){x}, function(y){0}, function(y){-1/6})
#     
#   print(solve(A, b) - answ(x[2:(length(x)-1)], y[2]))
  
  for (i in 2:(I-1)){
      A <- form.A.matrix.1(n, h1, I)
      b <- form.b.1(n, I, x, y[i], h1, f, g1, g2)
      
      res1[(i-1),] = TDMA.m(A, b)
  }
  
  for (j in 2:(J-1)){
    A <- form.A.matrix.1(n, h2, I)
    b <- form.b.1(m, I, y, x[j], h2, f, g3, g4)
    
    res2[,(j-1)] = TDMA.m(A, b)
  }
  
  image(x[2:(I-1)], y[2:(J-1)], t(res1+res2))
  title("experimental solve")
                          
}

form.A.matrix.1 <- function(n, h, I){
  A <- matrix(rep(0, times=n^2), n, n)
  
  d <- diag(A)
  diag.num <- -outer(seq(d),seq(d),"-")
  
  A[diag.num == 0] = rep(rep(-2/h^2), times=length(A[diag.num == 0]))
  A[diag.num == 1] = A[diag.num == -1] = rep(rep(1/h^2), times=length(A[diag.num == 1]))
  
  return (A)
}

form.b.1 <- function(n, I, x, y, h, f, g1, g2){
  b <- rep(0, times=n)
  
    for (i in 2:(I-1)){
      k = (i-2) + 1
      # !!! here we use f/2
      b[k] = -f(x[i], y)/2
      if (i == 2) b[k] = b[k] - g1(y)/h^2
      if (i == (I-1)) b[k] = b[k] - g2(y)/h^2
    }
  
  return (b)
}