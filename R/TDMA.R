# Script for testing Tridiagonal Matrix Algorythm

# @TODO: TDMA implementation 

#' @author Anton Shamin

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

TDMA.m <- function(A, d){
  b <- diag(A)
  diag.num <- -outer(seq(b),seq(b),"-")
  a = A[diag.num == -1]
  c = A[diag.num == 1]
  
  return (TDMA(a, b, c, d))
}