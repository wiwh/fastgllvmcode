#' Generates data from a GLLVM model.
#'
#' Returns a gllvm model of class fastglllvm with simulated data..
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
gen_fastgllvm <- function(n=100, p=5, q=1, k=0, family="gaussian", A=NULL, B=NULL, phi=NULL, X=NULL, Z=NULL, intercept=F){
  if(is.null(A)){
    if(is.null(p) | is.null(q)) stop("Either A, or p and q must be supplied.")
    A <- matrix(runif(p*q,-2, 2), p, q)
  } else {
    stopifnot(is.matrix(A))
    p <- nrow(A)
    q <- ncol(A)
  }
  if(is.null(B)){
    if(is.null(p) | is.null(k)) stop("Either B, or p and k must be supplied")
    if(k == 0){
      B <- matrix(0, p, 1)
    } else {
      B <- matrix(runif(p*k, -1, 1), p, k)  # this can be a matrix with 0 columns.
    }
  } else {
    stopifnot(is.matrix(B))
    stopifnot(nrow(B) == p)
    k <- ncol(B)
  }
  if(is.null(phi)){
    if(is.null(p)) stop("Either phi, or p must be supplied")
    phi <- rep(1, p)
  } else {
    stopifnot(length(phi) == p)
  }

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if(is.null(Z)) Z <- gen_Z(n, q)
  if(is.null(X)){
    if(k==0){
      X <- matrix(0, n, 1)
      k <- 1
    } else {
      X <- gen_X(n, k, intercept)
    }
  } else {
    if(intercept && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
  }
  if(intercept && k==0) warning("Intercept could not be added because k is set to 0. Set k=1 if intercept=True.")
  dat <- gen_Y(A, B, phi, Z, X, family)

  hist <- list(
    A = matrix(as.vector(A), 1, p*q),
    B = matrix(as.vector(B), 1, p*k),
    phi = matrix(phi, 1, p),
    crit = rep(0, 1)
  )

  fastgllvm <- structure(list(
    Y=dat$Y,
    natpar=dat$natpar,
    Z=Z,
    X=X,
    A=A,
    B=B,
    phi=phi,
    family=family,
    hist = hist,
    n=n,
    p=p,
    q=q,
    k=k
  ),class="fastgllvm")
  validate_fastgllvm(fastgllvm)
  fastgllvm
}

gen_Z <- function(n, q){
  matrix(stats::rnorm(n*q), n, q)
}


#' Returns a function used to generate Z. This is solely used within gllvm.
#'
#' @param method: one of "SA" (always random), "SP" (fixed)
#' @param H: number of draws to sample, usually set to 1 for SA and more for SP
#'
#' @return a function that returns a list of generated Z.
generate_Z_functionfactory <- function(n, q, H=1, method="SA"){
if(method == "SP"){
    Z_saved <- lapply(1:H, function(na) matrix(rnorm(n*q), n, q))
    generate_Z <- function(){
      Z_saved
    }
  }
  if(method == "SA"){
    generate_Z <- function(){
      lapply(1:H, function(na) matrix(rnorm(n*q), n, q))
    }
  }
  generate_Z
}


gen_X <- function(n, k, intercept){
  X <- matrix(rnorm(n*k), n, k)
  if(k >= 1 && intercept){
    X[,1] <- 1
  }
  X
}

gen_Y <- function(A, B, phi, Z, X, family){
  n <- nrow(Z)
  p <- nrow(A)
  q <- ncol(A)
  natpar <- compute_natpar(A, B, Z, X)
  Y <- sapply(1:p, function(j){
    switch(family$family,
      gaussian = rnorm(n, natpar[,j], phi[j]),
      binomial = rbinom(n, 1, family$linkinv(natpar[,j])),
      poisson  = rpois(n, family$linkinv(natpar[,j]))
    )
  })
  list(Y=Y, natpar=natpar)
}

compute_natpar <- function(A, B, Z, X){
  Z %*% t(A) + X %*% t(B)
}
