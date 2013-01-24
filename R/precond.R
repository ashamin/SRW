# TODO: 
#  1) Add SOR preconditioner
#  2) Rewrite minimal residuals method
###############################################################################

#' @description methods of solving differential schemes
#' @author Anton Shamin

gen.matrix.7diag.dominant <- function(n, up1, up2, lo1, lo2){
  
}

gen.matrix.diag.dominant <- function(n){
  A <- matrix(rnorm(n^2), n, n)
  diag(A) <- rnorm(n, mean=10)
  A
}

Jacobi.precond <- function(A){
  diag(diag(A))
}

Seidel.precond <- function(A){
  P <- A
  P[upper.tri(A)] = 0
  P
}

#Symmetric Successive Over Relaxation
# w from interval [.8, 1]
# best with w = 1
SSOR.precond.approach1 <- function(A, w){
  D <- diag(diag(A))
  L <- A
  U <- A
  L[upper.tri(L, diag=TRUE)]=0
  U[lower.tri(L, diag=TRUE)]=0
  (D/w + L) %*% ((w/(2-w))*solve(D)) %*% (D/w + U)
}

#Symmetric Successive Over Relaxation
# w from interval [.8, 1]
# best with w = 1
SSOR.precond <- function(A, w){
  D <- diag(diag(A))
  L <- A
  U <- A
  L[upper.tri(L, diag=TRUE)]=0
  U[lower.tri(L, diag=TRUE)]=0
  (D/w + L) %*% (solve(D)/w) %*% (D/w + U)
}

# w from interval [.8, 1]
# best with w = 1
relax.precond <- function(A, w){
  P <- A
  P[upper.tri(A)] = 0
  diag(P) <- diag(P)/w
  P
} 


slv <- function(A, f, P=Jacobi.precond(A), x0=rnorm(length(f)), epsilon=1e-5, maxit=200, tao=1, tao.compute=FALSE){
  it <- maxit
  repeat{
    it <-it - 1
    
    r <- f - A%*%x0
    delta <- solve(P)%*%r
    
    #minimal residuals method
    if (tao.compute){
      tao <- as.double(crossprod(A%*%r, r)/crossprod(A%*%r, A%*%r))
    }
    
    x <- x0 + tao*delta
    x0 <- x
    
    if (max(abs(f - A%*%x)) < epsilon) break 
    if (it == 0) stop("Iteration process obviously won't converge. \n Try to increase \"maxit\" value")
  }
  data.frame(root=x, residual=f - A%*%x, it=maxit-it)
}

test <- function(n, P=SSOR.precond.approach1(A, .8), x0=rnorm(length(f)), epsilon=1e-5, maxit=200, tao=1, tao.compute=FALSE){
  A <- gen.matrix.diag.dominant(n);
  f <- rnorm(n)
  slv(A, f, P)
}