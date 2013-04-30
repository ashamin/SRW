source("R/TDMA.R")

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
  diags <- form.diags(n);
  d <- rnorm(n) + 3;
  
  A <- form.matrix(diags$a, diags$b, diags$c)
  
  residuals1 <- solve(A, d) - TDMA.m(A, d)
  residuals2 <- solve(A, d) - TDMA(diags$a, diags$b, diags$c, d)
  View(residuals1)
  View(residuals2)
}