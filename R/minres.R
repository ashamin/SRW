
source("R/precond.R")
source("R/TDMA.R")

# Minimum residuals method. When tau.compute parameter equals to FALSE then function works as 
# simple iteration method with preconditioner P. Iteration parameter tau in this case is constant
# initial value of tau parameter (tau=1 as default)'
#
# MinRes method doesn't use preconditioner in iteration process. Preconditioned analog of MinRes is
# Minmum Correction method named here mincorr() - function.
minres <- function(A, f, P=Jacobi.precond(A), x0=rnorm(length(f)), epsilon=1e-5, maxit=200, tau=1, tau.compute=TRUE){
  it <- maxit
  step.residual <- double(0)
  iP = solve(P)
  repeat{
    it <-it - 1
    
    r <- (f - A%*%x0)
    
    # тао нулевое должно болтаться в районе 1!!!!!!!! и положительное! 1,5 2 может быть
    
    # if we compute tau we don't compute precoditioner and vice versa
    # this ussumption gave us not crusial number of iterations, but it's just 
    # experimental suggestion. not final solve of problem
    if (tau.compute) delta <- r
    else delta <- iP%*%r
    
    Ar = A%*%r
    
    #minimal residuals method
    if (tau.compute){
      tau <- as.double(crossprod(Ar, r)/crossprod(Ar, Ar))
      step.residual <- c(step.residual, norm(r))
    }
    
    x <- x0 + tau*delta
    x0 <- x
    
    if (max(abs(f - A%*%x)) < epsilon) break 
    if (it == 0) stop("Iteration process obviously won't converge. \n Try to increase \"maxit\" value")
  }
  return (list(root=x, residual=f - A%*%x, it=maxit-it, step.residual=step.residual))
}

# Mimum correction method. Preconditioned MinRes analog. Iterarion parameter tau compute differs 
# from MisRes method. 
mincorr <- function(A, f, P=Jacobi.precond(A), x0=rnorm(length(f)), epsilon=1e-5, maxit=200){
  
  #################################
#   taus = numeric(0)
  
  it <- maxit
  step.residual <- double(0)
  iP = solve(P)  
  repeat{
    it <-it - 1
    
    r <- (f - A%*%x0)
    w <- iP%*%r
    Aw = A%*%w
    
    tau <- as.double(crossprod(Aw, w)/crossprod(iP%*%Aw, Aw))
    step.residual <- c(step.residual, norm(r))
    
    #############################
#     taus = c(taus, tau)
    
    x <- x0 + tau*w
    x0 <- x
    
    if (max(abs(f - A%*%x)) < epsilon) break 
    if (it == 0) stop("Iteration process obviously won't converge. \n Try to increase \"maxit\" value")
  }
  ###############################
#   View(taus)
  return (list(root=x, residual=f - A%*%x, it=maxit-it, step.residual=step.residual))
}




# Parallel sequent experimental version of MINRES method
mincorr.par.seq <- function(A, f, P, x0=rnorm(length(f)), epsilon=1e-5, maxit=200){
  
  #################################
  #   taus = numeric(0)
  
  tau = 1
  
  it <- maxit
  step.residual <- double(0)
  
  iP = solve(P$Mx %*% P$My)
  
  n = sqrt(length(iP[1,]))

  repeat{
    it <-it - 1
    
    r <- (f - A%*%x0)    
    
    tmp = double(0)
    
    for (i in 1:n){
      k = ((i-1)*n + 1)
      p = ((i-1)*n + n)
      tmp = c(tmp, TDMA.m(P$Mx[k:p, k:p], r[k:p]))
    }
    
    corr = double(0)
    
    for (i in 1:n){
      k = ((i-1)*n + 1)
      p = ((i-1)*n + n)
      corr = c(corr, TDMA.m(P$My[k:p, k:p], tmp[k:p]))
    }
    
    Aw = A%*%corr
    
    tau <- as.double(crossprod(Aw, corr)/crossprod(iP%*%Aw, Aw))
    #tau = .95
    step.residual <- c(step.residual, norm(r))
    
    #############################
    #     taus = c(taus, tau)
    
    x <- x0 + tau*corr
    x0 <- x
    
    #print(max(abs(f - A%*%x)))
    #print(tau)
    
    if (max(abs(f - A%*%x)) < epsilon) break 
#     if (it == 0) stop("Iteration process obviously won't converge. \n Try to increase \"maxit\" value")
    if (it == 0) break;
  }
  ###############################
  #   View(taus)
  return (list(root=x, residual=f - A%*%x, it=maxit-it, step.residual=step.residual))
}


# Parallel sequent experimental version of MINRES method
mincorr.par <- function(A, f, P, x0=rnorm(length(f)), epsilon=1e-5, maxit=200){
  
  #################################
  #   taus = numeric(0)
  
  tau = 1
  
  it <- maxit
  step.residual <- double(0)
  
  iP = solve(P$Mx %*% P$My)
  
  n = sqrt(length(iP[1,]))
  
  repeat{
    it <-it - 1
    
    r <- (f - A%*%x0)    
    
    tmp = TDMA.m(P$Mx, r)
    
    tmp = double(0)
    
    for (i in 1:n){
      k = ((i-1)*n + 1)
      p = ((i-1)*n + n)
      tmp = c(tmp, TDMA.m(P$Mx[k:p, k:p], r[k:p]))
    }
    
    corr = double(0)
    
    for (i in 1:n){
      k = ((i-1)*n + 1)
      p = ((i-1)*n + n)
      corr = c(corr, TDMA.m(P$My[k:p, k:p], tmp[k:p]))
    }
    
    Aw = A%*%corr
    
    tau <- as.double(crossprod(Aw, corr)/crossprod(iP%*%Aw, Aw))
    #tau = .95
    step.residual <- c(step.residual, norm(r))
    
    #############################
    #     taus = c(taus, tau)
    
    x <- x0 + tau*corr
    x0 <- x
    
    #print(max(abs(f - A%*%x)))
    #print(tau)
    
    if (max(abs(f - A%*%x)) < epsilon) break 
    #     if (it == 0) stop("Iteration process obviously won't converge. \n Try to increase \"maxit\" value")
    if (it == 0) break;
  }
  ###############################
  #   View(taus)
  return (list(root=x, residual=f - A%*%x, it=maxit-it, step.residual=step.residual))
}



test <- function(n, P=SSOR.precond.approach1(A, 1), x0=rnorm(length(f)), epsilon=1e-5, maxit=500, tau=1, tau.compute=FALSE){
  A <- gen.matrix.diag.dominant(n);
  f <- rnorm(n)
  slv = minres(A, f, P, x0, epsilon, maxit, tau, tau.compute)
  print(slv$it)
}

minres.test <- function(n, P=SSOR.precond.approach1(A, 1), x0=rnorm(length(f)), epsilon=1e-5, maxit=10000){
  A <- gen.matrix.diag.dominant(n);
  f <- rnorm(n)
  slv = mincorr(A, f, P, x0, epsilon, maxit)
  #print(max(slv$residual))
  print(slv$it)
}