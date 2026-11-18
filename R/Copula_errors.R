#' @name Copula_errors
#' @title Copula Errors Function
#' @description This function calculates the Mean Squared Error (MSE) and relative error for Copula distributions.
#'
#' @usage Copula_errors(data, copula_type, param)
#'
#' @param data Matrix of data.
#' @param copula_type Character string specifying the type of Copula (e.g., "clayton", "gumbel", "frank", "t", "nested", "mixed", "galambos", "normal").
#' @param param Parameters for the specified Copula type.
#'
#' @return A named vector containing:
#' \item{MSE}{Mean Squared Error for the Copula distribution.}
#' \item{RelativeError}{Relative error for the Copula distribution.}
#'
#' @export
#'
#' @examples
#' # Example for Clayton Copula
#' if (!requireNamespace("copula", quietly = TRUE)) {
#'   utils::install.packages("copula")
#' }
#' library(copula)
#' set.seed(123)
#' n <- 1000
#' clayton.cop <- copula::claytonCopula(2, dim = 2)
#' data <- copula::rCopula(n, clayton.cop)
#' errors <- Copula_errors(data, "clayton", 2)
#' print(errors)
Copula_errors <- function(data, copula_type, param) {
  # Define Frobenius norm function
  frobenius.norm <- function(mat) {
    return(sqrt(sum(mat^2)))
  }

  # Generate Copula samples based on the type
  if (copula_type == "clayton") {
    cop <- copula::claytonCopula(param, dim = ncol(data))
  } else if (copula_type == "gumbel") {
    cop <- copula::gumbelCopula(param, dim = ncol(data))
  } else if (copula_type == "frank") {
    cop <- copula::frankCopula(param, dim = ncol(data))
  } else if (copula_type == "t") {
    cop <- copula::tCopula(param, dim = ncol(data), df = param)
  } else if (copula_type == "nested") {
    cop <- copula::onacopulaL("Clayton", list(2, 1, list(3, c(2, 3))))
  } else if (copula_type == "mixed") {
    cop <- copula::mixCopula(list(copula::gumbelCopula(2.5, dim = ncol(data)), copula::claytonCopula(3, dim = ncol(data))), w = c(0.4, 0.6))
  } else if (copula_type == "galambos") {
    cop <- copula::galambosCopula(param)
  } else if (copula_type == "normal") {
    rho_matrix <- matrix(c(1, param, param, 1), nrow = 2)
    cop <- copula::normalCopula(param = rho_matrix[lower.tri(rho_matrix)], dim = ncol(data), dispstr = "un")
  } else {
    stop("Unsupported Copula type")
  }

  # Generate samples from the Copula
  n <- nrow(data)
  U <- copula::rCopula(n, cop)

  # Calculate Mean Squared Error (MSE) and relative error
  MSE <- frobenius.norm(U - data)^2 / (n * ncol(data))
  RelativeError <- frobenius.norm(U - data)^2 / frobenius.norm(data)^2

  return(c(MSE = MSE, RelativeError = RelativeError))
}
