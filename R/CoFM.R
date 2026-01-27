#' @name CoFM
#' @title Generate Copula Factor Models Data and Perform Tests
#' @description This function simulates data based on a Copula Factor Model structure.
#' It generates factor scores, factor loadings, and error terms (using specified Copula
#' distributions), combines them to create the observed data, and then performs KMO and
#' Bartlett's tests on the generated data.
#'
#' @param n Integer. Sample size (number of rows). Default is 1000.
#' @param p Integer. Number of observed variables (columns). Default is 10.
#' @param m Integer. Number of factors. Default is 5.
#' @param type Character. The type of Copula for error terms. Options: "Clayton", "Gumbel", "Frank".
#' @param param Numeric. The parameter for the Copula (theta). Default is 2.0.
#'
#' @return A list containing:
#' \item{data}{The generated data matrix (n x p).}
#' \item{KMO}{The results of the Kaiser-Meyer-Olkin test.}
#' \item{Bartlett}{The results of Bartlett's test of sphericity.}
#' \item{True_Params}{A list containing the true parameters used (F, A, D, mu).}
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(123)
#'
#' # Clayton copula errors (toy size)
#' res1 <- CoFM::CoFM(n = 200, p = 6, m = 2, type = "Clayton", param = 2)
#' res1$KMO
#' res1$Bartlett
#'
#' # Gumbel copula errors (toy size)
#' res2 <- CoFM::CoFM(n = 150, p = 6, m = 2, type = "Gumbel", param = 2)
#' head(res2$data)
#'
#' @export
CoFM <- function(n = 1000, p = 10, m = 5, type = "Clayton", param = 2.0) {

  raw_mu <- stats::runif(p, 0, n)
  mu <- t(matrix(rep(raw_mu, n), p, n))

  mu0 <- as.matrix(stats::runif(m, 0, 1))
  sigma0 <- diag(stats::runif(m, 1, 2))
  F_scores <- matrix(MASS::mvrnorm(n, mu0, sigma0), nrow = n)

  A <- matrix(stats::runif(p * m, -1, 1), nrow = p)
  D <- as.matrix(diag(stats::runif(p, 0, 1)))

  cop_model <- switch(
    type,
    "Clayton" = copula::claytonCopula(param = param, dim = p),
    "Gumbel"  = copula::gumbelCopula(param = param, dim = p),
    "Frank"   = copula::frankCopula(param = param, dim = p),
    stop("Invalid Copula type. Please choose 'Clayton', 'Gumbel', or 'Frank'.")
  )

  U_sim <- copula::rCopula(n, cop_model)
  epsilon <- matrix(U_sim, nrow = n)

  data <- mu + F_scores %*% t(A) + epsilon

  cor_data <- stats::cor(data)
  kmo_res <- psych::KMO(data)
  bartlett_res <- psych::cortest.bartlett(cor_data, n = n)

  list(
    data = data,
    KMO = kmo_res,
    Bartlett = bartlett_res,
    True_Params = list(F = F_scores, A = A, D = D, mu = mu)
  )
}
