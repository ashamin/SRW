# TODO: 
#  1) Add SOR preconditioner
###############################################################################

#' @description preconditioners for solving linear systems
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