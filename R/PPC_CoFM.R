#' @name PPC.CoFM
#' @title PPC estimation for the Copula Factor Model
#' @description Projected-PCA (noise-reduction projection) for CoFM data.
#' @param data  n×p data matrix.
#' @param m     Number of factors.
#' @param A     True p×m loading matrix.
#' @param D     True p×p uniqueness matrix.
#' @return List with components
#'   \item{Ap2}{Estimated loadings.}
#'   \item{Dp2}{Estimated uniquenesses.}
#'   \item{MSESigmaA, MSESigmaD, LSigmaA, LSigmaD}{Error metrics.}
#' @examples
#' \donttest{
#'   X <- CoFM(n=1000, p=10, m=5, theta=2, distribution_type="Frank")
#'   res <- PPC.CoFM(X$data, m=5, X$A, diag(diag(X$D)))
#'   print(res)
#' }
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
PPC.CoFM <- function(data, m, A, D) {
  X <- scale(data)
  n <- nrow(X); p <- ncol(X)
  ## simple projection: remove first PC direction (empirical noise)
  P <- diag(n) - tcrossprod(rep(1/n, n))
  Xpro <- scale(P %*% X)
  Spro <- cov(Xpro)
  eig <- eigen(Spro)
  Ap2 <- eig$vectors[, 1:m] %*% diag(sqrt(eig$values[1:m]))
  Dp2 <- diag(Spro - diag(rowSums(Ap2^2)))

  MSESigmaA <- frobenius.norm(Ap2 - A)^2 / (p^2)
  MSESigmaD <- frobenius.norm(Dp2 - D)^2 / (p^2)
  LSigmaA   <- frobenius.norm(Ap2 - A)^2 / frobenius.norm(A)^2
  LSigmaD   <- frobenius.norm(Dp2 - D)^2 / frobenius.norm(D)^2
  list(Ap2 = Ap2, Dp2 = Dp2,
       MSESigmaA = MSESigmaA, MSESigmaD = MSESigmaD,
       LSigmaA   = LSigmaA,   LSigmaD   = LSigmaD)
}
