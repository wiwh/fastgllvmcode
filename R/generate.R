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
gen_gllvm <- function(n=100, p=NULL, q=NULL, k=0, family="normal", par=NULL, Z=NULL, X=NULL, scale=0.5){
  if(all(is.null(p), is.null(q), is.null(par)))
    stop("At least par, or both p and q, must be supplied.")
  if(!is.null(par)){
    p <- par$p
    q <- par$q
  } else {
    par <- gen_par(p, q, k, family, scale=scale)
  }


  if(is.null(Z)) Z <- gen_Z(n, q)
  if(is.null(X)) X <- gen_X(n, k)
  if(length(family)==1) family <- rep(family, p) # todo: at this point lead to a more efficient, dedicated function

  dat <- gen_Y(par, Z, X, family)

  list(Y=dat$Y, Z=Z, X=X, natpar=dat$natpar, par=par, n=n, p=p, q=q, k=k, scale=scale)
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
  natpar <- compute_natpar(par, Z, X)
  Y <- sapply(1:p, function(j){
    switch(family[j],
      normal = rnorm(n, natpar[,j], par$Psi[j]),
      bernoulli = rbinom(n, 1, sigmoid(natpar[,j])))
  })
  list(Y=Y, natpar=natpar)
}

#' Generates parameters for a gllvm
#'
#' returns a list of parameters.
#'
#'
gen_par <- function(p, q, k=0, family="normal", A=NULL, B=NULL, Psi=NULL, scale=0.5){
  # generate the (p, q) matrix of loadings
  if(is.null(A)) A <- matrix(rnorm(p*q), p, q)*scale
  # generate the vector of scale parameters
  if(is.null(Psi)) Psi <- runif(p, 0.5, 2)
  Psi[family=="bernoulli"] <- 1
  # generate the (p, k) matrix of coefficients
  if(is.null(B)){
    if(k==0){
      B <- matrix(0, p, 1)
    } else {
      B <- matrix(runif(p*k, -1, 1), p, k)
    }
  }
  list(A=A, B=B, Psi=Psi, family=family, p=p, q=q)
}

#' Modifies par with new parameters
#'
change_par <- function(par, A=NULL, B=NULL, Psi=NULL){
  if(!is.null(A)) par$A <- A
  if(!is.null(B)) par$B <- B
  if(!is.null(Psi)) par$Psi <- Psi
  par
}

compute_natpar <- function(par, Z, X){
  Z %*% t(par$A) + X %*% t(par$B)
}
