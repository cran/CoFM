#' The CoFM function is to generate Copula Factor Models data.
#'
#' The function supports various elliptical distributions for generating the data,
#' including: Elliptical-Normal Distribution, Elliptical-t Distribution, Elliptical Distribution.
#'
#' @name CoFM
#'
#' @usage CoFM(n, p, m, theta, distribution_type)
#'
#' @param n Sample size.
#' @param p Sample dimensionality.
#' @param m Number of factors.
#' @param theta A numerical parameter used in the copula distribution.
#' @param distribution_type The type of distribution.
#'
#' @return A list containing:
#' \item{data}{A matrix of generated data.}
#' \item{A}{A matrix representing the factor loadings.}
#' \item{D}{A diagonal matrix representing the unique variances.}
#'
#' @examples
#' library(copula)
#' library(MASS)
#' n <- 1000
#' p <- 10
#' m <- 5
#' theta <- 2.0
#' distribution_type <- "Clayton"
#' X <- CoFM(n, p, m, theta, distribution_type)
#'
#' @export
#'
#' @importFrom stats runif
#' @importFrom MASS mvrnorm
#' @importFrom copula claytonCopula gumbelCopula frankCopula rCopula
CoFM <- function(n, p, m, theta, distribution_type) {
  mu <- t(matrix(rep(runif(p, 0, 100), n), p, n))
  mu0 <- as.matrix(runif(m, 0))
  sigma0 <- diag(runif(m, 1))
  F <- matrix(MASS::mvrnorm(n, mu0, sigma0), nrow = n)
  A <- matrix(runif(p * m, -1, 1), nrow = p)
  D <- diag(runif(p, 1))

  if (distribution_type == "Clayton") {
    clayton.cop <- copula::claytonCopula(param = theta, dim = p)
    epsilon <- matrix(copula::rCopula(n, clayton.cop), nrow = n)
  } else if (distribution_type == "Gumbel") {
    gumbel.cop <- copula::gumbelCopula(param = theta, dim = p)
    epsilon <- matrix(copula::rCopula(n, gumbel.cop), nrow = n)
  } else if (distribution_type == "Frank") {
    frank.cop <- copula::frankCopula(param = theta, dim = p)
    epsilon <- matrix(copula::rCopula(n, frank.cop), nrow = n)
  } else {
    stop("Unsupported distribution_type")
  }

  data <- mu + F %*% t(A) + epsilon

  return(list(data = data, A = A, D = D))
}
