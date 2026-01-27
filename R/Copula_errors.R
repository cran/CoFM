#' @name Copula_errors
#' @title Generate Copula-Distributed Error Terms
#' @description
#' Generate random samples (error terms) from various copula distributions,
#' including Archimedean (Clayton, Gumbel, Frank), Elliptical (t, Normal),
#' Mixed, and Extreme-Value (Galambos) copulas. Useful for simulation studies
#' involving non-normal error structures.
#'
#' @param n Integer. The number of samples (rows) to generate.
#' @param type Character. The type of Copula to use. Options: "Clayton", "Gumbel",
#'   "Frank", "t", "Mixed", "Galambos", "Normal".
#' @param dim Integer. The dimension of the copula (number of columns/variables).
#'   Used for Archimedean/Elliptical copulas. For "Mixed" the default is 3.
#' @param param Numeric or Matrix. The main parameter for the copula (e.g., theta for
#'   Archimedean, correlation vector/matrix for Normal). If NULL, default values are used.
#' @param extra_params List. Additional parameters for specific copulas (e.g., `df` for t-Copula,
#'   `w` for Mixed).
#' @return A numeric matrix of dimension (n x dim) containing the generated random samples.
#'
#' @export
#'
#' @examples
#' # Examples should be fast and reproducible for CRAN checks
#' set.seed(123)
#'
#' # Example 1: Clayton Copula (toy example)
#' U_clayton <- Copula_errors(n = 200, type = "Clayton", dim = 2, param = 2)
#' head(U_clayton)
#'
#' # Example 2: t-Copula with degrees of freedom (toy example)
#' U_t <- Copula_errors(
#'   n = 200, type = "t", dim = 2, param = 0.7,
#'   extra_params = list(df = 4)
#' )
#' head(U_t)
#'
#' # Example 3: Multivariate Normal Copula (dim = 3)
#' # normalCopula() expects the upper-triangular correlations as a vector:
#' # (rho_12, rho_13, rho_23) for dim=3
#' rho_vec <- c(0.5, 0.3, 0.4)
#' U_normal <- Copula_errors(n = 200, type = "Normal", dim = 3, param = rho_vec)
#' head(U_normal)
Copula_errors <- function(n,
                          type = "Clayton",
                          dim = 2,
                          param = NULL,
                          extra_params = list()) {

  # Basic input checks (lightweight; avoids cryptic errors)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) {
    stop("'n' must be a single positive number.")
  }
  n <- as.integer(n)

  if (!is.numeric(dim) || length(dim) != 1L || is.na(dim) || dim <= 1) {
    stop("'dim' must be a single integer >= 2.")
  }
  dim <- as.integer(dim)

  if (!is.character(type) || length(type) != 1L) {
    stop("'type' must be a single character string.")
  }

  cop_model <- switch(
    type,

    "Clayton" = {
      p <- if (is.null(param)) 2 else param
      copula::claytonCopula(param = p, dim = dim)
    },

    "Gumbel" = {
      p <- if (is.null(param)) 3 else param
      copula::gumbelCopula(param = p, dim = dim)
    },

    "Frank" = {
      p <- if (is.null(param)) 5 else param
      copula::frankCopula(param = p, dim = dim)
    },

    "t" = {
      p <- if (is.null(param)) 0.7 else param
      dof <- if (is.null(extra_params$df)) 4 else extra_params$df
      copula::tCopula(param = p, dim = dim, df = dof)
    },

    "Mixed" = {
      # mixCopula needs copulas of the same dimension; default to dim = 3
      if (dim != 3L) {
        warning("For type = 'Mixed', dim is set to 3 to match the built-in mixture components.")
        dim <- 3L
      }
      g_cop <- copula::gumbelCopula(2.5, dim = 3)
      c_cop <- copula::claytonCopula(3, dim = 3)
      weights <- if (is.null(extra_params$w)) c(0.4, 0.6) else extra_params$w
      copula::mixCopula(list(g_cop, c_cop), w = weights)
    },

    "Galambos" = {
      # galambosCopula's dimension is defined by the object; set via dim argument
      p <- if (is.null(param)) 2 else param
      copula::galambosCopula(param = p, dim = dim)
    },

    "Normal" = {
      # normalCopula expects a VECTOR of correlations of length dim*(dim-1)/2 (dispstr = "un")
      if (is.null(param)) {
        # Default: equicorrelation 0.3 -> vector of required length
        rho <- 0.3
        param <- rep(rho, dim * (dim - 1L) / 2L)
      }

      # If user provided a correlation matrix, extract upper triangle
      if (is.matrix(param)) {
        if (nrow(param) != dim || ncol(param) != dim) {
          stop("If 'param' is a matrix, it must be a dim x dim correlation matrix.")
        }
        param <- param[upper.tri(param)]
      }

      needed <- dim * (dim - 1L) / 2L
      if (length(param) != needed) {
        stop(sprintf(
          "For type = 'Normal' with dim = %d, 'param' must have length %d (upper-triangular correlations).",
          dim, needed
        ))
      }

      copula::normalCopula(param = param, dim = dim, dispstr = "un")
    },

    stop("Invalid Copula type selected.")
  )

  copula::rCopula(n, cop_model)
}
