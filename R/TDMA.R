# Script for testing Tridiagonal Matrix Algorythm

# @TODO: TDMA implementation 

# Tridiagonal Matrix Algorythm
# a - lower diag, b - daig, c - upper diag
TDMA <- function(a, b, c, d){
  n <- length(d)
  x <- double(n);
  a <- c(0, a)
  c <- c(c, 0)
  
  c[1] <- c[1] / b[1]
  d[1] <- d[1] / b[1]
  
  for (i in 2:n){
    tmp <- b[i] - c[i-1]*a[i]
    c[i] <- c[i] / tmp
    d[i] <- (d[i] - d[i-1]*a[i]) / tmp
  }
  
  x[n] <- d[n]
  for (i in (n-1):1)
    x[i] <- d[i] - c[i]*x[i+1]
  
  return (x)
}

form.diags <- function(n){
  a <- rnorm(n-1) + 3
  b <- rnorm(n) + 10
  c <- rnorm(n-1) + 3
  
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
  diags <- form.diags(5);
  d <- rnorm(n) + 3;
  
  A <- form.matrix(diags$a, diags$b, diags$c)
  
  residuals <- solve(A, d) - TDMA(diags$a, diags$b, diags$c, d)
  View(residuals)
}