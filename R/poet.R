#' @name poet
#' @title POET: Principal Orthogonal complEment Thresholding
#' @description
#' Implements the POET method for large covariance matrix estimation (Fan, Liao & Mincheva, 2013).
#' The method assumes a factor model structure, estimates the low-rank component via PCA, and
#' applies thresholding to the sparse residual covariance matrix.
#'
#' @param X Numeric matrix (T x N). T is the number of time periods (rows), N is the number of variables (columns).
#' @param r Integer or NULL. User-specified number of factors. If NULL, r is estimated automatically using \code{method.r}.
#' @param r.max Integer. Upper bound for the number of factors when estimating r. Default is 10.
#' @param thresh Character. Thresholding type for the residual covariance. Options: "hard", "soft", "scad", "adapt".
#'   Default is "hard".
#' @param lambda Numeric or NULL. Thresholding parameter. If NULL, it defaults to \code{sqrt(log(N)/T)}.
#' @param gamma Numeric. Parameter for SCAD thresholding. Default is 3.7.
#' @param delta Numeric. Minimum eigenvalue bump to ensure positive definiteness of the residual covariance. Default is 1e-4.
#' @param method.r Character. Method to select the number of factors if r is NULL. Options: "IC1" (Bai & Ng, 2002)
#'   or "ER" (Eigenvalue Ratio, Ahn & Horenstein, 2013).
#' @return A list containing:
#' \item{Sigma.poet}{The estimated N x N POET covariance matrix.}
#' \item{Sigma.fact}{The estimated N x N low-rank (factor) covariance matrix.}
#' \item{Sigma.resid}{The estimated N x N thresholded residual covariance matrix.}
#' \item{F.hat}{Estimated factors (T x r).}
#' \item{Lambda.hat}{Estimated factor loadings (N x r).}
#' \item{r.hat}{The number of factors used (estimated or specified).}
#' \item{R.hat}{Same as Sigma.resid (for compatibility).}
#'
#' @references
#' Fan, J., Liao, Y., & Mincheva, M. (2013). Large covariance estimation by thresholding principal orthogonal complements.
#' \emph{Journal of the Royal Statistical Society: Series B}, 75(4), 603-680.
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(2025)
#' T_obs <- 40; N_var <- 15; r_true <- 2
#'
#' # Generate a simple factor model: X = F * Lambda' + U
#' Lambda <- matrix(stats::rnorm(N_var * r_true), N_var, r_true)
#' F_scores <- matrix(stats::rnorm(T_obs * r_true), T_obs, r_true)
#' U <- matrix(stats::rnorm(T_obs * N_var), T_obs, N_var)
#' X_sim <- F_scores %*% t(Lambda) + U  # T x N
#'
#' # Apply POET (choose r via IC1; use soft thresholding)
#' res <- poet(X_sim, r = NULL, method.r = "IC1", thresh = "soft")
#' res$r.hat
#' res$Sigma.poet[1:5, 1:5]
poet <- function(X, r = NULL, r.max = 10, thresh = "hard", lambda = NULL,
                 gamma = 3.7, delta = 1e-4, method.r = "IC1") {

  # ---------
  # Checks
  # ---------
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric T x N matrix.")
  T <- nrow(X)
  N <- ncol(X)
  if (T <= 1L || N <= 1L) stop("X must have at least 2 rows and 2 columns.")

  if (!is.null(r)) {
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r < 0) stop("'r' must be a single non-negative integer or NULL.")
    r <- as.integer(r)
  }

  if (!is.numeric(r.max) || length(r.max) != 1L || is.na(r.max) || r.max < 1) {
    stop("'r.max' must be a single integer >= 1.")
  }
  r.max <- as.integer(r.max)

  thresh <- match.arg(thresh, choices = c("hard", "soft", "scad", "adapt"))

  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
      stop("'lambda' must be a single non-negative number or NULL.")
    }
  }

  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || gamma <= 2) {
    stop("'gamma' must be a single number > 2 (SCAD parameter).")
  }
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) || delta < 0) {
    stop("'delta' must be a single non-negative number.")
  }

  if (!is.character(method.r) || length(method.r) != 1L || !(method.r %in% c("IC1", "ER"))) {
    stop("method.r must be 'IC1' or 'ER'.")
  }

  # Demean columns
  X <- scale(X, center = TRUE, scale = FALSE)


  if (is.null(r)) {
    r.max.use <- min(r.max, min(T, N) - 1L)
    if (r.max.use < 1L) {
      r <- 0L
    } else if (method.r == "IC1") {
      r <- select_r_IC1(X, r.max = r.max.use)
    } else {
      r <- select_r_ER(X, r.max = r.max.use)
    }
  }
  r.hat <- min(as.integer(r), min(T, N) - 1L)
  r.hat <- max(r.hat, 0L)

  if (r.hat == 0L) {
    # No factors: low-rank part is zero; residual is sample covariance of X
    Sigma.fact <- matrix(0, nrow = N, ncol = N)
    F.hat <- matrix(0, nrow = T, ncol = 0)
    Lambda.hat <- matrix(0, nrow = N, ncol = 0)
    E.hat <- X
  } else {
    svdX <- base::svd(X, nu = r.hat, nv = r.hat)
    U_r <- svdX$u[, 1:r.hat, drop = FALSE]
    V_r <- svdX$v[, 1:r.hat, drop = FALSE]
    D_r <- diag(svdX$d[1:r.hat], nrow = r.hat)

    # Factors and loadings (Fan et al. convention)
    F.hat <- sqrt(T) * U_r                       # T x r
    Lambda.hat <- V_r %*% D_r / sqrt(T)          # N x r

    Sigma.fact <- tcrossprod(Lambda.hat)         # N x N
    E.hat <- X - F.hat %*% t(Lambda.hat)         # T x N
  }

  S_resid <- crossprod(E.hat) / T                # N x N

  if (is.null(lambda)) {
    lambda <- sqrt(log(N) / T)
  }

  Sigma.resid <- threshold_matrix(S_resid, lambda = lambda, type = thresh, gamma = gamma)

  ev <- base::eigen(Sigma.resid, symmetric = TRUE)
  ev$values <- pmax(ev$values, delta)
  Sigma.resid <- ev$vectors %*% diag(ev$values, nrow = N) %*% t(ev$vectors)

  Sigma.poet <- Sigma.fact + Sigma.resid

  list(
    Sigma.poet = Sigma.poet,
    Sigma.fact = Sigma.fact,
    Sigma.resid = Sigma.resid,
    F.hat = F.hat,
    Lambda.hat = Lambda.hat,
    r.hat = r.hat,
    R.hat = Sigma.resid
  )
}

threshold_matrix <- function(M, lambda, type = "hard", gamma = 3.7) {
  M_sym <- (M + t(M)) / 2
  A <- M_sym
  n <- nrow(A)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      aij <- A[i, j]
      A[i, j] <- switch(
        type,
        "hard" = ifelse(abs(aij) > lambda, aij, 0),
        "soft" = sign(aij) * pmax(abs(aij) - lambda, 0),
        "scad" = scad_threshold(aij, lambda, gamma),
        "adapt" = {
          w <- ifelse(aij == 0, 1, 1 / abs(aij))
          ifelse(abs(aij) > lambda * w, aij, 0)
        }
      )
    }
  }
  (A + t(A)) / 2
}

scad_threshold <- function(x, lambda, a = 3.7) {
  absx <- abs(x)
  if (absx <= lambda) {
    0
  } else if (absx <= a * lambda) {
    sign(x) * (absx - lambda) / (1 - 1 / (a - 1))
  } else {
    x
  }
}

select_r_IC1 <- function(X, r.max = 10) {
  T <- nrow(X); N <- ncol(X)
  r.max <- min(r.max, min(T, N) - 1L)

  IC1 <- numeric(r.max)

  for (k in seq_len(r.max)) {
    svd_k <- base::svd(X, nu = k, nv = k)
    U_k <- svd_k$u[, 1:k, drop = FALSE]
    V_k <- svd_k$v[, 1:k, drop = FALSE]
    D_k <- diag(svd_k$d[1:k], nrow = k)

    F_k <- sqrt(T) * U_k
    Lambda_k <- V_k %*% D_k / sqrt(T)

    E_k <- X - F_k %*% t(Lambda_k)
    ssr <- sum(E_k^2) / (N * T)

    IC1[k] <- log(ssr) + k * (N + T) / (N * T) * log(N * T / (N + T))
  }

  which.min(IC1)
}

select_r_ER <- function(X, r.max = 10) {
  T <- nrow(X); N <- ncol(X)
  r.max <- min(r.max, min(T, N) - 1L)

  S <- crossprod(X) / T
  ev <- base::eigen(S, symmetric = TRUE)$values
  ev <- ev[1:(r.max + 1L)]
  ev <- pmax(ev, 1e-10)

  ratios <- ev[-length(ev)] / ev[-1]
  which.max(ratios)
}
