#' @name PPC_basic
#' @title Perform Basic Projected PCA (PPC) Estimation
#' @description
#' This function performs Projected Principal Component Analysis (PPC) to estimate factor loadings
#' and specific variances. It projects the data onto a specific subspace before performing eigen
#' decomposition. Unlike \code{\link{PPC.CoFM}}, this function does not calculate error metrics
#' against true parameters.
#'
#' @param data A matrix or data frame of input data (n x p).
#' @param m Integer. The number of principal components (factors) to extract.
#' @return A list containing:
#' \item{Apro}{Estimated projected factor loadings matrix (p x m).}
#' \item{Dpro}{Estimated projected uniquenesses vector (p).}
#' \item{Sigmahatpro}{The covariance matrix of the projected data.}
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(123)
#'
#' # 1. Generate toy data using CoFM
#' sim <- CoFM(n = 200, p = 6, m = 2, theta = 2.0, distribution_type = "Clayton")
#' obs_data <- sim$data
#'
#' # 2. Apply PPC method (extract 2 factors)
#' fit <- PPC_basic(data = obs_data, m = 2)
#'
#' # 3. Inspect estimates
#' head(fit$Apro)
#' fit$Dpro
PPC_basic <- function(data, m) {

  # 1. Input Validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m <= 0 || m > p) {
    stop("m must be a positive integer and cannot exceed the number of variables (ncol(data)).")
  }
  m <- as.integer(m)

  if (n < 2L || p < 2L) {
    stop("data must have at least 2 rows and 2 columns.")
  }

  # 2. Standardize Data
  X <- scale(X)

  w <- rep(c(0, 1), length.out = n)
  Xpro <- scale(X * w)  # row-wise scaling via recycling (each row multiplied by w[i])

  Sigmahatpro <- stats::cov(Xpro)

  eig <- base::eigen(Sigmahatpro, symmetric = TRUE)

  ind <- order(eig$values, decreasing = TRUE)
  lambdahat <- eig$values[ind]
  Q <- eig$vectors[, ind, drop = FALSE]

  Qhat <- Q[, 1:m, drop = FALSE]

  lam_use <- pmax(lambdahat[1:m], 0)
  Apro <- Qhat %*% diag(sqrt(lam_use), nrow = m)

  hpro <- diag(Apro %*% t(Apro))
  Dpro <- diag(Sigmahatpro) - hpro
  Dpro <- pmax(Dpro, 0)  # clamp small numerical negatives

  list(
    Apro = Apro,
    Dpro = Dpro,
    Sigmahatpro = Sigmahatpro
  )
}
