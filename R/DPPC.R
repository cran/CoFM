#' @name DPPC
#' @title Distributed Projection Principal Component Method
#' @description
#' Estimates a covariance matrix using a distributed projection principal
#' component method. Local loading matrices are aligned using singular
#' value decomposition before aggregation.
#'
#' @param data A numeric matrix or data frame with observations in rows
#'   and variables in columns.
#' @param m A positive integer specifying the number of principal components.
#' @param n1 An integer specifying the sample size of each local subset.
#' @param K A positive integer specifying the number of local subsets.
#'
#' @return Abarpro, Dbarpro, Sigmahatbarpro.
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10; m <- 3
#' F <- matrix(rnorm(n * m), nrow = n, ncol = m)
#' A <- matrix(runif(p * m, -1, 1), nrow = p, ncol = m)
#' E <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
#' X <- F %*% t(A) + E
#'
#' fit <- DPPC(data = X, m = 3, n1 = 50, K = 2)
#' dim(fit$Abarpro)
#' dim(fit$Dbarpro)
#' dim(fit$Sigmahatbarpro)
#'
#' @importFrom stats cov
#' @export
DPPC <- function(data, m, n1, K) {

  X <- as.matrix(data)
  n <- nrow(X); p <- ncol(X); m <- as.integer(m)
  stopifnot(m >= 1, m <= p, n1 >= 2, n1 <= n, K >= 1)

  X1pro <- scale(X)
  Apro <- list(); Dpro <- list(); Sigmahatpro <- list()
  for (i in seq_len(K)) {
    L <- sample.int(n, size = n1, replace = FALSE)
    Xi <- X1pro[L, , drop = FALSE]
    Xpro <- scale(Xi, center = TRUE, scale = TRUE)
    stopifnot(all(is.finite(Xpro)))

    Sigmahatproi <- cov(Xpro)
    Sigmahatproi <- (Sigmahatproi + t(Sigmahatproi)) / 2

    eigi <- base::eigen(Sigmahatproi, symmetric = TRUE)
    ind <- order(eigi$values, decreasing = TRUE)
    lambdahat <- pmax(eigi$values[ind][1:m], 0)
    Q <- eigi$vectors[, ind, drop = FALSE]
    Qhat <- Q[, 1:m, drop = FALSE]

    Aproi <- Qhat %*% diag(sqrt(lambdahat), nrow = m)
    hproi <- diag(Aproi %*% t(Aproi))
    dij <- pmax(diag(Sigmahatproi) - hproi, 0)
    Dhatproi <- diag(dij, nrow = p)

    Apro[[i]] <- Aproi
    Dpro[[i]] <- Dhatproi
    Sigmahatpro[[i]] <- Sigmahatproi
  }
  A.ref <- Apro[[1]]
  for (k in seq_len(K)) {
    sv <- svd(t(Apro[[k]]) %*% A.ref)
    H <- sv$u %*% t(sv$v)
    Apro[[k]] <- Apro[[k]] %*% H
  }
  Abarpro <- Reduce("+", Apro) / K
  Dbarpro <- Reduce("+", Dpro) / K
  Sigmahatbarpro <- Reduce("+", Sigmahatpro) / K
  Sigmahatbarpro <- (Sigmahatbarpro + t(Sigmahatbarpro)) / 2

  return(list(Abarpro = Abarpro, Dbarpro = Dbarpro,
              Sigmahatbarpro = Sigmahatbarpro))
}
