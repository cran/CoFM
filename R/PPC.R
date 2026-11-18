#' @name PPC
#' @title Projection-based Perturbation Principal Component Analysis
#' @param data is a total  data set.
#' @param m is the number of principal component
#' @return Apro, Dpro, Sigmahatpro
#' @export
#' @examples
#' library(copula)
#' library(MASS)
#' n=1000
#' p=10
#' m=5
#' mu=t(matrix(rep(runif(p,0,n),n),p,n))
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#' F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
#' A=matrix(runif(p*m,-1,1),nrow=p)
#' D=as.matrix(diag(rep(runif(p,0,1))))
#' theta = 2.0
#' clayton.cop <- claytonCopula(param = theta, dim = p)
#' U_clayton <- rCopula(n, clayton.cop)
#' epsilon_clayton <- qnorm(U_clayton)
#' data=mu+F%*%t(A)+epsilon_clayton
#' PPC(data=data,m=5)
#' @export
PPC <- function(data, m) {
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)
  m <- as.integer(m)
  Xpro <- scale(X)
  Sigmahatpro <- cov(Xpro)
  eig <- eigen(Sigmahatpro, symmetric = TRUE)
  ind <- order(eig$values, decreasing = TRUE)
  lambdahat <- pmax(eig$values[ind][1:m], 0)
  Q <- eig$vectors[, ind, drop = FALSE]
  Qhat <- Q[, 1:m, drop = FALSE]
  Apro <- Qhat %*% diag(sqrt(lambdahat), nrow = m)
  hpro <- diag(Apro %*% t(Apro))
  Dpro <- diag(Sigmahatpro) - hpro
  Dpro <- pmax(Dpro, 0)
  Dpro <- diag(Dpro, nrow = p)
  list(Apro = Apro,Dpro = Dpro,Sigmahatpro = Sigmahatpro)
}
