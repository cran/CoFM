#' @name PC.CoFM
#' @title PC estimation for the Copula Factor Model
#' @description Classical PCA solution for CoFM data.
#' @param data  n×p data matrix.
#' @param m     Number of factors.
#' @param A     True p×m loading matrix.
#' @param D     True p×p uniqueness matrix.
#' @return List with components
#'   \item{A2}{Estimated loadings.}
#'   \item{D2}{Estimated uniquenesses.}
#'   \item{MSESigmaA, MSESigmaD, LSigmaA, LSigmaD}{Error metrics.}
#' @examples
#' \donttest{
#'   X <- CoFM(n=1000, p=10, m=5, theta=2, distribution_type="Gumbel")
#'   res <- PC.CoFM(X$data, m=5, X$A, diag(diag(X$D)))
#'   print(res)
#' }
#' @export
PC.CoFM <- function(data, m, A, D) {
  X  <- scale(data)
  R  <- cor(X)
  p  <- nrow(R)
  eig <- eigen(R)
  A2 <- eig$vectors[, 1:m] %*% diag(sqrt(eig$values[1:m]))
  D2 <- diag(R - diag(rowSums(A2^2)))

  MSESigmaA <- frobenius.norm(A2 - A)^2 / (p^2)
  MSESigmaD <- frobenius.norm(D2 - D)^2 / (p^2)
  LSigmaA   <- frobenius.norm(A2 - A)^2 / frobenius.norm(A)^2
  LSigmaD   <- frobenius.norm(D2 - D)^2 / frobenius.norm(D)^2
  list(A2 = A2, D2 = D2,
       MSESigmaA = MSESigmaA, MSESigmaD = MSESigmaD,
       LSigmaA   = LSigmaA,   LSigmaD   = LSigmaD)
}
