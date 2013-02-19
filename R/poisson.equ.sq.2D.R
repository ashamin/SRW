# implementation of solving Dirihle promlem for Poisson Equation in square area
#
#   delta u = d^2u(x, y)/dx^2 + d^2u(x, y)/dy^2 = -f(x, y)
#
# Boundary conditions:
#
#   u(0, y) = g1(y), u(a, y) = g2(y), y from [0, b]
#
#   u(x, 0) = g3(x), u(x, b) = g4(x), x from [0, a]
#
# g1..g4 must agree with conditions
#
#   g1(0) = g3(0)
#   g1(b) = g4(0)
#   g2(0) = g3(a)
#   g2(b) = g4(a)
#
#' @author Anton Shamin


# TODO: implement case when I != J

init <- function(){
  source("R/precond.R")
}

# Test area
# You may change this for solving your own Poisson equation by default
#

#Test 1
p <- 19
a.test <- 1
b.test <- 1

g1.test <- function(y) -2*y - 4*y^2
g2.test <- function(y) 4 - 12*y -4*y^2
g3.test <- function(x) x + 3*x^2
g4.test <- function(x) 3*x^2 - 9*x - 6

f.test <- function(x, y)  2 + 2*pi^2*(p - 15)*sin(pi*x)*sin(pi*y)

right.answer <- function(x, y) x + 3*x^2 - 10*x*y - 2*y - 4*y^2 + (p-15)*sin(pi*x)*sin(pi*y)
#
#
#



slv.poisson <- function(a=a.test, b=b.test, g1=g1.test, g2=g2.test,
                g3=g3.test, g4=g4.test, f=f.test,
                I=20, J=20, h1=a/(I-1), h2=b/(J-1)){
  
  x <- seq(from=0, to=a, by=h1)
  y <- seq(from=0, to=b, by=h2)
  
  n <- (I-2)*(J-2)
  
  A <- form.A.matrix(n, h1, h2, I, J)
    
  b <- form.b(n, I, J, x, y, h1, h2, f, g1, g2, g3, g4)
  
#   ra <- matrix(rep(0, length(x)*length(y)), length(x), length(y), byrow=TRUE)
#   for (i in 1:length(x))
#     for (j in 1:length(y))
#       ra[i, j] = right.answer(x[i], y[j])
#   image(x, y, ra)
  
  res = data.frame(solve = slv(A, b, P=SSOR.precond.approach1(A, .8), maxit=1000)$root, 
                  answer = get.res.vector(n, I, J, x, y))
  
  View(res)
  
  print(h1)
  print(h2)
  print(as.double(A[1,]%*%res$answer))
  print(b[1])
  View(A)
  View(b)
}

form.A.matrix <- function(n, h1, h2, I, J){
  A <- matrix(rep(0, times=n^2), n, n)
  
  d <- diag(A)
  diag.num <- -outer(seq(d),seq(d),"-")
  A[diag.num == 0] = rep(-(2/h1^2 + 2/h2^2), times=length(A[diag.num == 0]))
  A[diag.num == 1] = A[diag.num == -1] = rep(1/h1^2, times=length(A[diag.num == 1]))
  A[diag.num == (I-2)] = A[diag.num == -(J-2)] = rep(1/h2^2, times=length(A[diag.num == (I-2)]))
  
  d <- A[diag.num == 1]
  d[seq(from=(I-2), to=length(d), by=(I-2))] = 0
  A[diag.num == 1] = A[diag.num == -1] = d
  
  return (A)
}

form.b <- function(n, I, J, x, y, h1, h2, f, g1, g2, g3, g4){
  b <- rep(0, times=n)
  
  for (j in 2:(J-1))
    for (i in 2:(I-1)){
      k = (I-2)*(j-2) + (i-2) + 1
      b[k] = -f(x[i], y[j])
      if (j == 2) b[k] = b[k] - g3(x[i])/h2^2
      if (j == (J-1)) b[k] = b[k] - g4(x[i])/h2^2
      if (i == 2) b[k] = b[k] - g1(y[j])/h1^2
      if (i == (I-1)) b[k] = b[k] - g2(y[j])/h1^2
    }
  
  return (b)
}

get.res.vector <- function(n, I, J, x, y){
  res <- rep(0, times=n)
  
  for (j in 2:(J-1))
    for (i in 2:(I-1)){
      k = (I-2)*(j-2) + (i-2) + 1
      res[k] = right.answer(x[i], y[j])
    }
  
  return (res)
}