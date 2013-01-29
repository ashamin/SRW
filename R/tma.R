# Script for testing Tridiagonal Matrix Algorythm

#solves matrix with 3 diagonals with sweep method
#@parameters a - lower diag, b - diag, c - upper diag, d - vector
sweep <- function(a, b, c, d){
  n <- length(d)
  x <- real(n)
  alpha <- c(rep(0, times=n+1))
  betta <- c(rep(0, times=n+1))
  a <- c(0, a)
  c <- c(c, 0)
  alpha[2:(n+1)] <- c[1:n] / (b[1:n] - a[1:n]*alpha[1:n])
  betta[2:(n+1)] <- (a[1:n]*betta[1:n]-d[1:n]) / (b[1:n]-a[1:n]*alpha[1:n])
  
  x[n] <- betta[n+1]
  
  x[(n-1):1] <- alpha[n:2]*x[n:2] - betta[n:2]
  x
}

TDMA.solve <- function(a, b, c, d){
  N <- length(d)
  alpha <- real(N+1)
  betta <- real(N+1)
  x <- real(N)
  
  a = c(0, a)
  c = c(c, 0)
  
  alpha[2] = -c[1] / b[1]
  betta[2] = d[1] / b[1]
  
  for (i in 2:N){
    tmp = a[i]*alpha[i] + b[i]
    alpha[i] = -c[i] / tmp
    betta[i] = (d[i] - a[i]*betta[i]) / tmp
  }
  
  #View(alpha)
  
  #x[N] = betta[N]
  x[N] = (d[N] - a[N]*betta[N]) / (b[N] + a[N]*alpha[N])
  for (i in (N-1):1)
    x[i] = alpha[i+1]*x[i+1] + betta[i+1]
  
  return (x)
}

TDMA.wiki <- function(a, b, c, d){
  N <- length(d)
  alpha <- double(N)
  betta <- double(N)
  x <- double(N)

  a <- c(0, a)
  c <- c(c, 0)
  
  c[1] <- c[1] / b[1]
  d[1] <- d[1] / b[1]
  
  for (n in 2:N){
    temp <- b[n] - a[n]
    c[n] <- c[n] / temp
    d[n] <- (d[n] - a[n]*d[n-1]) / temp
  }
  
  x[N] <- d[N]
  for (n in (N-1):1)
    x[n] <- d[n] - c[n] * x[n+1]
            
  return (x)  
}

form.diags <- function(n){
  a <- rnorm(n-1)
  b <- rnorm(n)
  c <- rnorm(n-1)
  
  return (list(a=a, b=b, c=c))
}

form.matrix <- function(a, b, c){
  M <- diag(b)
  diag.num <- -outer(seq(b),seq(b),"-")
  M[diag.num == -1] = a
  M[diag.num == 1] = c
  
  return (M)
}

test <- function(n){
  diags <- form.diags(n)
  diags$b <- diags$b + 10
  diags$a <- diags$a + 3
  diags$c <- diags$c + 3
  
  A <- form.matrix(diags$a, diags$b, diags$c)
  d <- rnorm(n)
  
  #residuals <- solve(A, d) - TDMA.solve(diags$a, diags$b, diags$c, d)
  View(A)
  View(solve(A, d))
  residuals.gauss <- A%*%solve(A, d) - d
  View(residuals.gauss)
  TDMA.result <- TDMA.wiki(diags$a, diags$b, diags$c, d)
  View(TDMA.result)
  residuals.TDMA.solve <- A%*%TDMA.result - d
  View(residuals.TDMA.solve)
  #View(residuals)
}

