#' @name POET
#' @title Projection-on-Orthogonal-Complement Estimation and Transformation
#' @param X A numeric n by p data matrix, with observations in rows
#'   and variables in columns.
#' @param r.max A positive integer specifying the maximum number of
#'   factors considered for selection. The effective upper bound is
#'   limited by the sample size and number of variables.
#' @param delta A non-negative numeric value added to the diagonal
#'   of the residual covariance estimate. Default is 0.01.
#' @param eps A positive numeric tolerance used to stabilize
#'   eigenvalue-based calculations and residual variances.
#'   Default is 1e-8.
#' @return Ahat, Dhat, Sigmahat, Fhat.
#' @references
#' Fan, J., Liao, Y. and Mincheva, M. (2013).
#' Large covariance estimation by thresholding principal orthogonal complements.
#' \emph{Journal of the Royal Statistical Society: Series B},
#' 75(4), 603--680.
#'
#' @examples
#' set.seed(123)
#' n <- 100; p <- 10; r <- 3
#' F <- matrix(rnorm(n * r), nrow = n, ncol = r)
#' A <- matrix(runif(p * r, -1, 1), nrow = p, ncol = r)
#' E <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
#' X <- F %*% t(A) + E
#' res <- POET(X = X, r.max = 5, delta = 0.01, eps = 1e-8)
#' dim(res$Sigmahat)
#' res$rhat
#' Sigma.true <- A %*% t(A) + diag(0.25, p)
#' MSESigma <- mean((res$Sigmahat - Sigma.true)^2)
#' MSESigma
#' @export
POET <- function(X, r.max = 10, delta = 0.01, eps = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (!is.numeric(X) || any(!is.finite(X))) stop("X must be a finite numeric matrix.")
  if (n < 2 || p < 2) stop("X must have at least 2 rows and 2 columns.")
  if (length(r.max) != 1 || !is.finite(r.max) || r.max < 1 || r.max != floor(r.max)) stop("r.max must be a positive integer.")

  Xc <- scale(X, center = TRUE, scale = FALSE)
  r.max <- min(as.integer(r.max), n - 1, p - 1)
  Sigma.c <- cov(Xc)
  Sigma.c <- (Sigma.c + t(Sigma.c)) / 2
  eig <- eigen(Sigma.c, symmetric = TRUE)
  ind <- order(eig$values, decreasing = TRUE)
  vals <- pmax(eig$values[ind], 0)
  Q <- eig$vectors[, ind, drop = FALSE]
  E.k <- sum(Xc^2) - c(0, cumsum(colSums((Xc %*% Q)^2)))[seq_len(r.max + 1)]
  V.k <- pmax(E.k / (n * p), .Machine$double.eps)
  penalty.k <- (0:r.max) * (n + p) * log(n * p / (n + p)) / (n * p)
  IC1.values <- log(V.k) + penalty.k
  rhat <- which.min(IC1.values[-1])
  Q.r <- Q[, seq_len(rhat), drop = FALSE]
  lambda.r <- vals[seq_len(rhat)]
  Lambda.sqrt <- diag(sqrt(lambda.r), nrow = rhat)
  Ahat <- Q.r %*% Lambda.sqrt
  Fhat <- Xc %*% Q.r %*% diag(1 / sqrt(pmax(lambda.r, eps)), nrow = rhat)
  Sigma.factor <- Ahat %*% t(Ahat)
  Ehat <- Xc - (Xc %*% Q.r) %*% t(Q.r)
  Rhat.raw <- crossprod(Ehat) / (n - 1)
  tau <- sqrt(log(p) / n)
  Dhat <- Rhat.raw * (abs(Rhat.raw) >= tau)
  diag(Dhat) <- pmax(diag(Rhat.raw), eps) + delta
  Sigmahat <- Sigma.factor + Dhat
  Sigmahat <- (Sigmahat + t(Sigmahat)) / 2
  list(Ahat = Ahat, Dhat = Dhat, Sigmahat = Sigmahat, Fhat = Fhat, rhat = rhat)
}
