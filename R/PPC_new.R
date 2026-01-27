#' @name PPC_new
#' @title Center-then-PCA: Projection on the Orthogonal Complement of the Mean Vector
#' @description
#' This function performs a specific type of Projected PCA where the data is projected onto
#' the orthogonal complement of the mean vector. It effectively applies the centering projection
#' \eqn{P = I - (1/n)J} (where \eqn{J} is the all-ones matrix), optionally rescales the columns,
#' and then performs PCA on the covariance matrix. This allows estimation of factor loadings and
#' residual variances after removing the mean structure.
#'
#' @param data A matrix or data frame of input data (n x p).
#' @param m Integer. The number of principal components (factors) to keep.
#'
#' @return A list containing:
#' \item{Apro}{Estimated factor loading matrix (p x m).}
#' \item{Dpro}{Estimated residual variances (p x p diagonal matrix).}
#' \item{Sigmahatpro}{Covariance matrix of the projected data.}
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(1)
#' dat <- matrix(stats::rnorm(200), ncol = 4)
#' ans <- PPC_new(data = dat, m = 2)
#' str(ans)
#' head(ans$Apro)
PPC_new <- function(data, m) {

  # 1. Input Validation
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
    stop("m must be a positive integer and cannot exceed the number of variables.")
  }
  m <- as.integer(m)

  # 2. Projection onto orthogonal complement of the mean vector
  # P %*% X where P = I - (1/n)J is exactly column-centering:
  # X_centered = X - 1 * colMeans(X)
  Xpro <- scale(X, center = TRUE, scale = TRUE)

  # 3. Covariance of the projected data
  Sigmahatpro <- stats::cov(Xpro)

  # 4. Eigen decomposition (covariance matrix is symmetric)
  eig <- base::eigen(Sigmahatpro, symmetric = TRUE)

  # Sort descending (defensive)
  ind <- order(eig$values, decreasing = TRUE)
  lambdahat <- eig$values[ind]
  Q <- eig$vectors[, ind, drop = FALSE]

  # Extract top m eigenvectors
  Qhat <- Q[, 1:m, drop = FALSE]

  # 5. Estimate loadings: Apro = Q_m * sqrt(lambda_m)
  lam_use <- pmax(lambdahat[1:m], 0)
  Apro <- Qhat %*% diag(sqrt(lam_use), nrow = m)

  # 6. Estimate uniquenesses (diagonal residual variances)
  hpro <- diag(Apro %*% t(Apro))
  Dpro_vec <- diag(Sigmahatpro) - hpro
  Dpro_vec <- pmax(Dpro_vec, 0)  # clamp tiny negatives from numerical error
  Dpro <- diag(Dpro_vec, nrow = p)

  # 7. Return
  list(
    Apro = Apro,
    Dpro = Dpro,
    Sigmahatpro = Sigmahatpro
  )
}
