# implementation of solving boundary problem for heat equation
#
# du(t, x)/dt = a^2 * d^2u(t, x)/dx^2 + f(x, t)
#
# 0 < x < l, 0 < t < T,
#
# u(0, x) = u0(x)
#
# Boundary conditions:
#
# alpha1 * du(t, 0)/dx = beta1 * u(t, 0) - mu1(t)
#
# -alpha2 * du(t, l)/dx = beta2 * u(t, l) - mu2(t)
#
#' @author Anton Shamin

source("R/TDMA.R")
library("animation")


# Test area
# You may change this for solving your own heat equation by default
#

#Test 1
# y.limit = 2
# p <- 19
# z <- function(t) p*t/16
# 
# a.test <- sqrt(p)/4
# l.test <- pi/2
# Tr.test <- 1
# 
# u0.test <- function(x) sin(x)
# 
# f.test <- function(t, x) p/16*exp(z(t))*sin(x)
# 
# alpha1.test <- 0
# beta1.test <- 1
# mu1.test <- function(t) 0
# 
# alpha2.test <- 1
# beta2.test <- 0
# mu2.test <-function(t) 0
# 
# right.answer <- function(t, x) ((exp(z(t)) + exp(-z(t)))/2)*sin(x)

#Test 2
# y.limit = 1
# p <- 3
# z <- function(t) p*t/2
# 
# a.test <- sqrt(p/2)
# l.test <- pi/2
# Tr.test <- 1
# 
# u0.test <- function(x) cos(x)
# 
# f.test <- function(t, x) p/2*exp(-z(t))*cos(x)
# 
# alpha1.test <- 3
# beta1.test <- 1
# mu1.test <- function(t) (1+z(t))*exp(-z(t))
# 
# alpha2.test <- 0
# beta2.test <- 1
# mu2.test <-function(t) 0
# 
# right.answer <- function(t, x) (1+z(t))*exp(-z(t))*cos(x)


#Test 3
# y.limit = 1
# p <- 8
# 
# a.test <- 1
# l.test <- 1
# Tr.test <- 1
# 
# u0.test <- function(x) sin(sqrt(p)*x)
# 
# f.test <- function(t, x) 0
# 
# alpha1.test <- 1
# beta1.test <- 1
# mu1.test <- function(t) -sqrt(p)*exp(-p*t)
# 
# alpha2.test <- 0
# beta2.test <- 1
# mu2.test <-function(t) exp(-p*t)*sin(sqrt(p))
# 
# right.answer <- function(t, x) exp(-p*t)*sin(sqrt(p)*x)

#Test 4
# y.limit = 10
# p <- 13
# s <- (p+1)/p
# 
# a.test <- 1
# l.test <- 1
# Tr.test <- 1
# 
# u0.test <- function(x) x*exp(sqrt(s)*x)
# 
# f.test <- function(t, x) -2*sqrt(s)*exp(sqrt(s)*x + s*t)
# 
# alpha1.test <- 0
# beta1.test <- 1
# mu1.test <- function(t) 0
# 
# alpha2.test <- -1
# beta2.test <- 1 + sqrt(s)
# mu2.test <-function(t) 0
# 
# right.answer <- function(t, x) x*exp(sqrt(s)*x+s*t)

#Test 5
y.limit = .5
p <- 23
q <- p^(1/6)

a.test <- 1
l.test <- 1
Tr.test <- 1

u0.test <- function(x) x*cos(q*x)

f.test <- function(t, x) 2*q*exp(-q^2*t)*sin(q*x)

alpha1.test <- -1
beta1.test <- 0
mu1.test <- function(t) exp(-q^2*t)

alpha2.test <- 0
beta2.test <- 1
mu2.test <-function(t) exp(-q^2*t)*cos(q)

right.answer <- function(t, x) x*exp(-q^2*t)*cos(q*x)
#
#
#







#' @TODO!! what about N & K. may be troubles with indices!
slv <- function(Tr=Tr.test, l=l.test, a=a.test, u0=u0.test, f=f.test,
                alpha1=alpha1.test, beta1=beta1.test, mu1=mu1.test,
                alpha2=alpha2.test, beta2=beta2.test, mu2=mu2.test,
                N=50, K=50, sigma=1, tau=Tr/(N-1), h=l/(K-1),
                p.mode="spectrogram", view.values=FALSE){
  
  x <- seq(from=0, to=l, by=h)
  t <- seq(from=0, to=Tr, by=tau)
  y <- matrix(rep(0, times=length(t)*length(x)), length(t), length(x), byrow=TRUE)
  
  if (sigma > 1 | sigma < 0) stop("sigma value must lies in [0, 1]")
  else{
    if (sigma != 0 & sigma != 1){
      mu1.m <- function(n) mu1(t[n]+tau/2) + (alpha1*h / (2*a^2)) * f(t[n]+tau/2, 0)
      mu2.m <- function(n) mu2(t[n]+tau/2) + (alpha2*h / (2*a^2)) * f(t[n]+tau/2, l)
    }
    else{
      mu1.m <- function(n) mu1(t[n+1]) + (alpha1*h / (2*a^2)) * f(t[n+1], 0)
      mu2.m <- function(n) mu2(t[n+1]) + (alpha2*h / (2*a^2)) * f(t[n+1], l)
    }
  }
  
  phi <- function(n, k) f(t[n] + sigma*tau, x[k])
  
  #defining helpful functions corresponds to alpha1 value
  if (alpha1 != 0){
    delta1 = sigma*(alpha1 + h*beta1) + (alpha1*h^2)/(2*tau*a^2)
    xi1 = alpha1*sigma/delta1
    
    nu1 <- function(n) (1/delta1) * (h*mu1.m(n) + (1-sigma)*alpha1*y[n, 2] +
                                       y[n, 1]*(alpha1*h^2 / (2*tau*a^2) - (1-sigma)*(alpha1+beta1*h)))
  }
  else{
    xi1 = 0
    nu1 <- function(n) mu1(t[n+1])/beta1
  }
  
  #defining helpful functions corresponds to alpha2 value
  if (alpha2 != 0){
    delta2 = sigma*(alpha2 + h*beta2) + (alpha2*h^2)/(2*tau*a^2)
    xi2 = alpha2*sigma/delta2
    
    nu2 <- function(n) (1/delta2) * (h*mu2.m(n) + (1-sigma)*alpha2*y[n, K-1] +
                                       y[n, K]*(alpha2*h^2 / (2*tau*a^2) - (1-sigma)*(alpha2+beta2*h)))
  }
  else{
    xi2 = 0
    nu2 <- function(n) mu2(t[n+1])/beta2
  }
  
  A = C = sigma*a^2
  B = 2*sigma*a^2 + (h^2)/tau
  
  D <- function(n, k) (2*(1-sigma)*a^2 - h^2/tau)*y[n, k] - (1-sigma)*(y[n, k-1] + y[n, k+1])*a^2 -
    phi(n, k)*h^2
  
  #filling first layer with u0 fuction values
  y[1, 1:K] = u0(x[1:K])
  
  D.s <- double(K-2)
  
  #compute all layers
  for (n in 1:(N-1)){
    
    #bulding scheme
    for (i in 1:(K-2))
      D.s[i] <- D(n, i+1)
    
    scheme <- scheme.build(xi1=xi1, nu1=nu1(n), xi2=xi2, nu2=nu2(n), A, B, C, D.s, y[n+1,])
    y[n+1,] = TDMA(a=scheme$a, b=scheme$b, c=scheme$c, d=scheme$d)
  }
  
  ra <- matrix(rep(0, length(t)*length(x)), length(t), length(x), byrow=TRUE)
  for (i in 1:length(t))
    for (j in 1:length(x))
      ra[i, j] = right.answer(t[i], x[j])
  
  
  if (p.mode == "spectrogram"){
    
    par(mfcol=c(1, 2))
    image(t, x, y)
    title("solve")
    image(t, x, ra)
    title("right answer")
  }
  else if (p.mode == "animation"){
    
    opts <- ani.options(interval=.2, nmax=N+1)
    
    for (i in 1:N){
      plot(x, y[i,], type="l", col="blue", xlim=c(0, l), ylim=c(0, y.limit))
      lines(x, ra[i,], type="l", col="magenta")
      
      ani.pause()
    }
    
    ani.pause()
  }
  
  print("Error")
  print(norm(ra - y) / norm(y))
#   print("r")
#   print(a^2 * tau / h)

  if (view.values == TRUE){
    View(y)
    View(ra)
  }
}

scheme.build <- function(xi1, nu1, xi2, nu2, A, B, C, D, y.n.plus1){
  K <- length(y.n.plus1)
  a <- double(K-1); b <- double(K); c <- double(K-1); d <- double(K)
  
  a[K-1] = -xi2; c[1] = -xi1;
  b[1] = 1; b[K] = 1
  d[1] = nu1; d[K] = nu2
  
  a[1:(K-2)] = A; c[2:(K-1)] = C
  b[2:(K-1)] = -B
  d[2:(K-1)] = D[1:(K-2)]
  
  return (list(a=a, b=b, c=c, d=d))
}


