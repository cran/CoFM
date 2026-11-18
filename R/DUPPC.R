#' @name DUPPC
#' @title Distributed U-projection Principal Component Method
#'
#' @description
#' Estimates a covariance matrix using a distributed U-projection principal
#' component method. The data are projected onto the orthogonal complement
#' of a direction vector, and local factor loading matrices are aligned
#' using singular value decomposition before aggregation.
#'
#' @param data A numeric n by p matrix or data frame, with observations
#'   in rows and variables in columns. All values must be finite.
#' @param m A positive integer specifying the number of principal components.
#' @param n1 An integer specifying the number of observations sampled for
#'   each local subset, between 2 and the total sample size.
#' @param K A positive integer specifying the number of local subsets.
#' @param seed An optional integer random seed for reproducibility.
#'   Default is NULL, which leaves the random number generator unchanged.
#' @param u An optional nonzero projection direction vector of length n.
#'   If NULL, the first principal component score vector is used.
#'
#' @return A list containing: Abarpro, Dbarpro, Sigmahatbarpro.
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10; m <- 3
#' F <- matrix(rnorm(n * m), nrow = n, ncol = m)
#' A <- matrix(runif(p * m, -1, 1), nrow = p, ncol = m)
#' E <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
#' X <- F %*% t(A) + E
#' result1 <- DUPPC(data = X, m = 3, n1 = 50, K = 2, seed = 123)
#' dim(result1$Abarpro)
#' dim(result1$Dbarpro)
#' dim(result1$Sigmahatbarpro)
#' u <- F[, 1]
#' result2 <- DUPPC(data = X, m = 3, n1 = 50, K = 2, seed = 123, u = u)
#' dim(result2$Sigmahatbarpro)
#'
#' @importFrom stats cov sd
#' @export
DUPPC <- function(data, m, n1, K, seed = NULL, u = NULL) {
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)
  m <- as.integer(m)
  X0 <- scale(X)
  eig <- eigen(cov(X0), symmetric = TRUE)

  if (is.null(u)) {
    u <- as.vector(X0 %*% eig$vectors[, 1])
  }

  u <- as.vector(u)
  nu <- sqrt(sum(u^2))
  u <- u / nu
  P <- diag(n) - tcrossprod(u)
  X1pro <- scale(P %*% X0)
  Apro <- list(); Dpro <- list(); Sigmahatpro <- list()

  for (i in 1:K) {
    idx <- sample.int(n, n1, replace = FALSE)
    R <- matrix(0, n1, n); R[cbind(1:n1, idx)] <- 1
    X1 <- R %*% X1pro
    X1_mean <- colMeans(X1); X1_sd <- apply(X1, 2, sd)
    X1_sd[!is.finite(X1_sd) | X1_sd == 0] <- 1
    Xpro <- scale(X1, center = X1_mean, scale = X1_sd)
    Sigmahatproi <- cov(Xpro); Sigmahatproi <- (Sigmahatproi + t(Sigmahatproi)) / 2
    Sigmahatpro[[i]] <- Sigmahatproi
    eig1 <- eigen(Sigmahatproi, symmetric = TRUE)
    lambdahat <- pmax(eig1$values[1:m], 0)
    Q <- eig1$vectors[, 1:m, drop = FALSE]
    Aproi <- Q %*% diag(sqrt(lambdahat), m, m)
    hproi <- diag(Aproi %*% t(Aproi))
    dij <- pmax(diag(Sigmahatproi) - hproi, 0)
    Dhatproi <- diag(dij, nrow = p)
    Apro[[i]] <- Aproi; Dpro[[i]] <- Dhatproi
  }

  A.ref <- Apro[[1]]
  for (k in 1:K) { sv <- svd(t(Apro[[k]]) %*% A.ref); H <- sv$u %*% t(sv$v); Apro[[k]] <- Apro[[k]] %*% H }

  Abarpro <- Reduce("+", Apro) / K
  Dbarpro <- Reduce("+", Dpro) / K
  Sigmahatbarpro <- Reduce("+", Sigmahatpro) / K
  Sigmahatbarpro <- (Sigmahatbarpro + t(Sigmahatbarpro)) / 2
  list(Abarpro = Abarpro, Dbarpro = Dbarpro, Sigmahatbarpro = Sigmahatbarpro)
}
