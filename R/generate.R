#' Generates data from a GLLVM model.
#'
#' Returns data of type gllvm.
#'
#' @param n the number of observations to draw
#' @param p the number of manifest variables
#' @param q the number of latent variables
#' @param k the number of covariates common to all responses, k=0 yields an intercept of 0, while k=1 yields a random intercept for each response, and k>1 generates covariates.
#' @param family one of "normal", "bernoulli", or a vector of both of size p, specifying individual responses' distribution
#' @param par optional parameters in a list
#' @param sigma sigma
#'
#' @return a list corresponding to the model
#' @export
gen_gllvm <- function(n, p, q, k=0, family="normal", par=NULL, Z=NULL, X=NULL){
  if(is.null(Z)) Z <- gen_Z(n, q)
  if(is.null(X)) X <- gen_X(n, k)
  if(length(family)==1) family <- rep(family, p)
  if(is.null(par)) par <- gen_par(p, q, k, family)
  Y <- gen_Y(par, Z, X, family)

  list(Y=Y, Z=Z, X=X, par=par)
}

gen_Z <- function(n, q){
  matrix(stats::rnorm(n*q), n, q)
}

gen_X <- function(n, k){
  if(k<=1){
    X <- matrix(1, n, 1)
  } else {
    X <- matrix(stats::rnorm(n*k), n, k)
    X[,1] <- 1
  }
  X
}

gen_Y <- function(par, Z, X, family){
  n <- nrow(Z)
  p <- nrow(par$A)
  natpar <- Z %*% t(par$A) + X %*% t(par$B)
  Y <- sapply(1:p, function(j){
    switch(family[j],
      normal = rnorm(n, natpar[,j], par$psi[j]),
      bernoulli = rbinom(n, 1, sigmoid(natpar[,j])))
  })
  Y
}

#' Generates parameters for a gllvm
#'
#' returns a list of parameters
#'
#'
gen_par <- function(p, q, k=0, family){
  # generate the (p, q) matrix of loadings
  A <- matrix(rnorm(p*q), p, q)
  # generate the vector of scale parameters
  psi <- runif(p, 0.5, 2)
  psi[family=="bernoulli"] <- 1
  # generate the (p, k) matrix of coefficients
  if(k==0){
    B <- matrix(0, p, 1)
  } else {
    B <- matrix(rnorm(p*k), p, k)
  }
  list(A=A, Psi=psi, B=B)
}


compute_natpar <- function(par, Z, X){
  Z %*% t(par$A) + X %*% t(par$B)
}
