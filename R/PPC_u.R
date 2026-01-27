#' @name PPC_u
#' @title Projection-on-Complement PCA (Generalized)
#' @description
#' Projects the data onto the orthogonal complement of a given vector \code{u}, eliminating the
#' effect of \code{u}, and then performs PCA on the projected data. This is useful for removing
#' specific trends (e.g., time trends, common market factors) before analysis.
#'
#' @param data A matrix or data frame of input data (n x p).
#' @param m Integer. Number of principal components to retain.
#' @param u Numeric vector of length n. The projection direction to be removed from the data.
#'   Will be normalized internally.
#'
#' @return A list containing:
#' \item{Apro}{Estimated factor loading matrix (p x m).}
#' \item{Dpro}{Estimated residual variances (p x p diagonal matrix).}
#' \item{Sigmahatpro}{Covariance matrix of the projected data.}
#' \item{u}{The normalized projection vector used.}
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(123)
#' dat <- matrix(stats::rnorm(200), ncol = 4)
#' u0 <- seq_len(nrow(dat))     # e.g., a linear trend to remove
#' res <- PPC_u(data = dat, m = 2, u = u0)
#' res$u
#' head(res$Apro)
PPC_u <- function(data, m, u) {

  # 1. Input checks
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  if (n < 2L || p < 2L) {
    stop("data must have at least 2 rows and 2 columns.")
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m <= 0 || m > p) {
    stop("m must be a positive integer and cannot exceed the number of variables (ncol(data)).")
  }
  m <- as.integer(m)

  if (missing(u)) stop("projection vector u must be supplied.")
  u <- as.numeric(u)
  if (length(u) != n) stop("length of u must equal number of rows in data.")
  if (anyNA(u)) stop("u must not contain NA values.")

  nu <- sqrt(sum(u^2))
  if (!is.finite(nu) || nu == 0) stop("u must be a non-zero finite vector.")
  u <- u / nu

  Xs <- scale(X)

  proj_part <- tcrossprod(u, drop(crossprod(u, Xs)))
  Xpro <- scale(Xs - proj_part)

  Sigmahatpro <- stats::cov(Xpro)

  eig <- base::eigen(Sigmahatpro, symmetric = TRUE)
  ind <- order(eig$values, decreasing = TRUE)
  lambdahat <- eig$values[ind]
  Q <- eig$vectors[, ind, drop = FALSE]

  Qhat <- Q[, 1:m, drop = FALSE]

  lam_use <- pmax(lambdahat[1:m], 0)
  Apro <- Qhat %*% diag(sqrt(lam_use), nrow = m)

  hpro <- diag(Apro %*% t(Apro))
  Dpro_vec <- diag(Sigmahatpro) - hpro
  Dpro_vec <- pmax(Dpro_vec, 0)           # clamp tiny negatives
  Dpro <- diag(Dpro_vec, nrow = p)

  list(
    Apro = Apro,
    Dpro = Dpro,
    Sigmahatpro = Sigmahatpro,
    u = u
  )
}
