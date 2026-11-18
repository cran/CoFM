#' @name UPPC
#' @title Projection-on-Complement Principal Component Analysis
#' @description
#' Performs principal component analysis after projecting the
#' standardized data onto the orthogonal complement of the
#' leading principal component score direction.
#'
#' @param data A numeric n by p data matrix, with observations
#'   in rows and variables in columns.
#' @param m A positive integer specifying the number of
#'   principal components to retain. Must be smaller than p.
#' @return Apro, Dpro, Sigmahatpro.
#'
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10; m <- 3
#' F <- matrix(rnorm(n * m), nrow = n)
#' A <- matrix(runif(p * m, -1, 1), nrow = p)
#' E <- matrix(rnorm(n * p, sd = 0.5), nrow = n)
#' X <- F %*% t(A) + E
#' result <- UPPC(data = X, m = 3)
#' dim(result$Apro)
#' length(result$Dpro)
#' dim(result$Sigmahatpro)
#' @export
UPPC <- function(data, m) {
  X <- as.matrix(data)
  n <- nrow(X); p <- ncol(X)
  m <- as.integer(m)
  X0 <- scale(X)
  if (any(!is.finite(X0)))
    stop("data contains constant columns; standardization failed.")
  Sigma0 <- cov(X0)
  eig0 <- eigen(Sigma0, symmetric = TRUE)
  u <- as.vector(X0 %*% eig0$vectors[, 1])
  nu <- sqrt(sum(u^2))
  u <- u / nu
  P <- diag(n) - tcrossprod(u)
  Xproj <- P %*% X0
  Xpro <- scale(Xproj)
  if (any(!is.finite(Xpro)))
    stop("Projected data contains constant columns.")
  Sigmahatpro <- cov(Xpro)
  Sigmahatpro <- (Sigmahatpro + t(Sigmahatpro)) / 2
  eig <- eigen(Sigmahatpro, symmetric = TRUE)
  ind <- order(eig$values, decreasing = TRUE)
  lambdahat <- pmax(eig$values[ind], 0)
  Q <- eig$vectors[, ind, drop = FALSE]
  Qhat <- Q[, seq_len(m), drop = FALSE]
  lam <- lambdahat[seq_len(m)]
  Apro <- Qhat %*% diag(sqrt(lam), nrow = m)
  hpro <- rowSums(Apro^2)
  Dpro <- pmax(diag(Sigmahatpro) - hpro, 0)
  list(
    Apro = Apro,
    Dpro = Dpro,
    Sigmahatpro = Sigmahatpro
  )
}
