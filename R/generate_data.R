# Random Numbers Generation -----------
# TODO:  write these in c++!

gen_norm <- function(linpar, phi, nobs, p){
  t(matrix(rnorm(p*nobs), p, nobs)*phi) + linpar
}

gen_poisson <- function(linpar, nobs, p){
  sapply(1:p, function(j)rpois(nobs, exp(linpar[,j])))
}
gen_binom <- function(linpar, size, nobs, p){
  sapply(1:p, function(j) rbinom(nobs, size[j], 1/(1+exp(-linpar[,j]))))
}

gen_binom_fast <- function(linpar, size, nobs, p){ # GEN BINOM FAST!
  if(all(size==1)){
    (matrix(runif(nobs*p), nobs, p) < 1/(1+exp(-linpar)))*1
  } else {
  tot <- sum(size)
  runifs <- matrix(runif(tot*nobs), nobs, tot)
  probs <- 1/(1+exp(-linpar))
  size.cumsum <- c(1,cumsum(size))
  sapply(1:p, function(j){
    if(size.cumsum[j+1] - size.cumsum[j] == 1){
      (runifs[,size.cumsum[j]] < probs[j])*1
    } else {
      rowSums(runifs[, size.cumsum[j]:size.cumsum[j+1]]<probs[j]) # TODO: change to colSums, which is faster for R!?
    }
  })
  }
}

gen_Z <- function(nobs, q){
  matrix(stats::rnorm(nobs*q), nobs, q)
}

#' Returns a function used to generate Z. This is solely used within gllvm.
#'
#' @param method: one of "SA" (always random), "SP" (fixed)
#' @param H: number of draws to sample, usually set to 1 for SA and more for SP
#'
#' @return a function that returns a list of generated Z.
generate_Z_functionfactory <- function(nobs, q, H=1, method="SA"){
if(method == "SP"){
    Z_saved <- lapply(1:H, function(na) matrix(rnorm(nobs*q), nobs, q))
    generate_Z <- function(){
      Z_saved
    }
  }
  if(method == "SA"){
    generate_Z <- function(){
      lapply(1:H, function(na) matrix(rnorm(nobs*q), nobs, q))
    }
  }
  generate_Z
}

gen_X <- function(nobs, k, intercept){
  X <- matrix(rnorm(nobs*k), nobs, k)
  if(k >= 1 && intercept){
    X[,1] <- 1
  }
  X
}

#' Generates responses from a gllvm, given all parameters and values
#'
#' Returns a matrix Y
#' @param A: loadings
#' @param B: beta
#' @param phi: scale parameters
#' @param Z: nobs x q matrix of latent variables
#' @param X: nobs x k design matrix
#' @param families: a named list with names, c("gaussian", "binomial", "poisson", "negbin"), of indices corresponding to the family of the responses.
#' @param known_par: a named list with names, c("gaussian", "binomial", "poisson", "negbin"), of known parameters (used for the size "n" of binomial responses), appearing in the same order as the corresponding element in `families`
#' @return a list corresponding to the model
generate_y <- function(linpar, phi, families, nobs, p, q) {
  Y <- matrix(NA,nobs,p)
  if(length(families$id$gaussian)>0){
    Y[,families$id$gaussian] <- gen_norm(linpar$linpar[,families$id$gaussian], phi=phi[families$id$gaussian], nobs=nobs, p=length(families$id$gaussian))
  }
  if(length(families$id$binomial)>0){
    Y[,families$id$binomial] <- gen_binom(linpar$linpar[,families$id$binomial], size=rep(1,p), nobs=nobs, p=length(families$id$binomial))
  }
  if(length(families$id$poisson)>0){
    Y[,families$id$poisson] <- gen_poisson(linpar$linpar[,families$id$poisson], nobs=nobs, p=length(families$id$poisson))
  }
  Y
}
