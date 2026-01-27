#' @name FanPC_basic
#' @title Perform Basic FanPC Factor Analysis
#' @description
#' This function performs factor analysis using a principal-component (FanPC) approach.
#' It estimates the factor loading matrix and uniquenesses from the correlation matrix
#' of the input data. Unlike \code{\link{FanPC_CoFM}}, this function does not calculate
#' error metrics against true parameters, making it suitable for simple estimation tasks.
#'
#' @param data A matrix or data frame of input data (n x p).
#' @param m Integer. The number of principal components (factors) to extract.
#' @return A list containing:
#' \item{AF}{Estimated factor loadings matrix (p x m).}
#' \item{DF}{Estimated uniquenesses vector (p).}
#' \item{SigmahatF}{The correlation matrix of the input data.}
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(123)
#'
#' # 1. Generate synthetic data using CoFM (toy example)
#' sim <- CoFM(n = 200, p = 6, m = 2, type = "Clayton", param = 2.0)
#' obs_data <- sim$data
#'
#' # 2. Apply FanPC method (extract 2 factors)
#' fit <- FanPC_basic(data = obs_data, m = 2)
#'
#' # 3. Inspect estimates
#' head(fit$AF)  # Estimated loadings
#' fit$DF        # Estimated uniquenesses
FanPC_basic <- function(data, m) {

  # 1. Input Validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }

  X <- as.matrix(data)
  p <- ncol(X)

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m <= 0 || m > p) {
    stop("m must be a positive integer and cannot exceed the number of variables (ncol(data)).")
  }
  m <- as.integer(m)

  # 2. Standardize Data
  X <- scale(X)

  # 3. Correlation Matrix
  SigmahatF <- stats::cor(X)

  # 4. Eigen Decomposition
  eig <- base::eigen(SigmahatF, symmetric = TRUE)

  # Sort eigenvalues/vectors descending (defensive)
  idx <- order(eig$values, decreasing = TRUE)
  lambdahat <- eig$values[idx]
  Q <- eig$vectors[, idx, drop = FALSE]

  # 5. Extract Loadings (PCA-style): loadings = eigenvectors * sqrt(eigenvalues)
  # Ensure non-negative eigenvalues for numerical stability
  lam_use <- pmax(lambdahat[1:m], 0)
  AF <- Q[, 1:m, drop = FALSE] %*% diag(sqrt(lam_use), nrow = m)

  # 6. Estimate Uniquenesses
  # Communalities = diag(AF %*% t(AF))
  hF <- diag(AF %*% t(AF))

  # Since SigmahatF is a correlation matrix, diag(SigmahatF) = 1
  DF <- diag(SigmahatF) - hF

  # Small negative values may occur due to numeric issues; clamp at 0
  DF <- pmax(DF, 0)

  # 7. Return Results
  list(
    AF = AF,
    DF = DF,
    SigmahatF = SigmahatF
  )
}
