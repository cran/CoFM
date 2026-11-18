
#' @name DPOET
#' @title Distributed POET covariance matrix estimation
#'
#' @description
#' Estimates a high-dimensional covariance matrix using a distributed Principal
#' Orthogonal complEment Thresholding (POET) procedure. Local subsets are randomly
#' sampled from the full dataset, and their covariance estimates are averaged.
#'
#' @param data A numeric matrix or data frame. Rows represent observations and
#'   columns represent variables. All values must be finite.
#' @param r.max A positive integer specifying the maximum number of latent factors
#'   considered by the information criterion. The effective maximum is limited to
#'   min(r.max, n1 - 1, p - 1), where p is the number of variables. Default is 10.
#' @param n1 An integer between 2 and the total sample size. The number of
#'   observations randomly sampled without replacement for each local subset.
#' @param K A positive integer specifying the number of local subsets. Each subset
#'   is sampled independently from the full dataset.
#' @param delta A non-negative finite numeric value added to the diagonal entries
#'   of the thresholded residual covariance matrix. Default is 0.01.
#' @param eps A positive finite numeric value used as a numerical lower bound for
#'   factor eigenvalues and residual covariance diagonal entries. Default is 1e-8.
#'
#' @return A list containing:
#' \describe{
#'   \item{Abar}{A p by r_common matrix obtained by averaging the zero-padded
#'     local factor loading matrices, where r_common is the maximum local factor number.}
#'   \item{Dbar}{A p by p matrix representing the average thresholded residual covariance estimate.}
#'   \item{Sigmahatbar}{A p by p symmetric matrix representing the average local POET covariance estimate.}
#'   \item{rbar}{The mean of the estimated local factor numbers.}
#' }
#'
#' @details
#' For each local subset, the number of latent factors is selected using an
#' information criterion. The residual covariance matrix is estimated using
#' hard thresholding with threshold sqrt(log(p) / n1).
#' Local subsets may overlap because they are independently sampled.
#' For reproducibility, call set.seed() before running DPOET().
#'
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10; r <- 3
#' F <- matrix(rnorm(n * r), nrow = n, ncol = r)
#' A <- matrix(runif(p * r, -1, 1), nrow = p, ncol = r)
#' U <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
#' X <- F %*% t(A) + U
#' result1 <- DPOET(data = X, r.max = 3, n1 = 50, K = 2, delta = 0.01, eps = 1e-8)
#' dim(result1$Sigmahatbar)
#' result1$rbar
#' result2 <- DPOET(data = X, r.max = 5, n1 = 50, K = 2)
#' dim(result2$Abar)
#' result2$rbar
#' @importFrom stats cov rnorm runif
#' @export
DPOET <- function(data, r.max = 10, n1, K, delta = 0.01, eps = 1e-8) {
  if (!is.matrix(data) && !is.data.frame(data)) stop
  X <- as.matrix(data)
  if (!is.numeric(X) || any(!is.finite(X))) stop

  n <- nrow(X); p <- ncol(X)
  if (n < 2 || p < 2) stop
  valid_integer <- function(x, lower, upper = Inf) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x <= upper && x == floor(x)
  }

  if (!valid_integer(r.max, 1)) stop
  if (!valid_integer(K, 1, .Machine$integer.max)) stop
  if (!valid_integer(n1, 2, n)) stop
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) stop
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) stop
  r.max.local <- as.integer(min(r.max, n1 - 1, p - 1))
  K <- as.integer(K); n1 <- as.integer(n1)
  Apro <- vector("list", K); Dpro <- vector("list", K)
  Sigmahatpro <- vector("list", K); rpro <- integer(K)
  for (i in seq_len(K)) {

    idx <- sample(seq_len(n), n1, replace = FALSE)
    X1 <- X[idx, , drop = FALSE]
    Xc <- scale(X1, center = TRUE, scale = FALSE)
    Sigma.c <- cov(Xc)
    Sigma.c <- (Sigma.c + t(Sigma.c)) / 2
    eig <- eigen(Sigma.c, symmetric = TRUE)
    ind <- order(eig$values, decreasing = TRUE)
    vals <- pmax(eig$values[ind], 0)
    Q <- eig$vectors[, ind, drop = FALSE]
    projected <- Xc %*% Q[, seq_len(r.max.local), drop = FALSE]
    E.k <- sum(Xc^2) - c(0, cumsum(colSums(projected^2)))
    V.k <- pmax(E.k / (n1 * p), .Machine$double.eps)
    penalty.k <- (0:r.max.local) * (n1 + p) * log(n1 * p / (n1 + p)) / (n1 * p)
    IC1.values <- log(V.k) + penalty.k
    rhat <- which.min(IC1.values[-1L])
    Q.r <- Q[, seq_len(rhat), drop = FALSE]
    lambda.r <- vals[seq_len(rhat)]
    Lambda.sqrt <- diag(sqrt(lambda.r), nrow = rhat)
    Ahat <- Q.r %*% Lambda.sqrt
    Sigma.factor <- Ahat %*% t(Ahat)
    Ehat <- Xc - (Xc %*% Q.r) %*% t(Q.r)
    Rhat.raw <- crossprod(Ehat) / (n1 - 1)
    tau <- sqrt(log(p) / n1)
    Dhat <- Rhat.raw * (abs(Rhat.raw) >= tau)
    diag(Dhat) <- pmax(diag(Rhat.raw), eps) + delta
    Sigmahat <- Sigma.factor + Dhat
    Sigmahat <- (Sigmahat + t(Sigmahat)) / 2
    Apro[[i]] <- Ahat; Dpro[[i]] <- Dhat
    Sigmahatpro[[i]] <- Sigmahat; rpro[i] <- rhat
  }
  r_common <- max(rpro)
  Apro <- lapply(Apro, function(A) {
    if (ncol(A) < r_common) {
      cbind(A, matrix(0, nrow = nrow(A), ncol = r_common - ncol(A)))
    } else {
      A
    }
  })
  Abar <- Reduce("+", Apro) / K
  Dbar <- Reduce("+", Dpro) / K
  Sigmahatbar <- Reduce("+", Sigmahatpro) / K
  Sigmahatbar <- (Sigmahatbar + t(Sigmahatbar)) / 2
  list(Abar = Abar, Dbar = Dbar, Sigmahatbar = Sigmahatbar, rbar = mean(rpro))
}
