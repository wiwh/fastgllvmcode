#' Generate multivariate Normal data conditional on the linear parameter
gen_norm <- function(linpar, phi, nobs, p){
  if(any(phi < 0)) stop("Negative values for the variance parameters.")
  t(matrix(rnorm(p*nobs), p, nobs)*sqrt(phi)) + linpar
}

#' Generate multivariate Poisson data conditional on the linear parameter
gen_poisson <- function(linpar, nobs, p){
  sapply(1:p, function(j)rpois(nobs, exp(linpar[,j])))
}

#' Generate multivariate binomial data conditional on the linear parameter
gen_binom <- function(linpar, size, nobs, p){
  sapply(1:p, function(j) rbinom(nobs, size[j], 1/(1+exp(-linpar[,j]))))
}

#' Generate multivariate binomial data conditional on the linear parameter
#'
#' This is experimental and is faster than `gen_binom`.
gen_binom_fast <- function(linpar, size, nobs, p){
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


#' Generates responses from a gllvm.
#'
#' Returns a matrix Y
#' @inheritParams gllvmprime
#' @return a list containing the responses and the various objects generated.

gen_Y <- function(Z, X, parameters, families, linpar=NULL, XB=NULL) {
  # compute linpar based on the provided information
  if (is.null(linpar)) {
    if(is.null(Z)) {
      Z <- gen_Z(nrow(X), ncol(parameters$A))
    }
    linpar <- compute_linpar(Z, parameters$A, X, parameters$B, XB)
  }
  Y <- matrix(NA,nrow(linpar$linpar), ncol(linpar$linpar))

  if(length(families$id$gaussian)>0){
    Y[,families$id$gaussian] <- gen_norm(linpar$linpar[,families$id$gaussian,drop=F], phi=parameters$phi[families$id$gaussian], n=nrow(linpar$linpar), p=length(families$id$gaussian))
  }
  if(length(families$id$binomial)>0){
    Y[,families$id$binomial] <- gen_binom(linpar$linpar[,families$id$binomial,drop=F], size=rep(1,length(families$id$binomial)), n=nrow(linpar$linpar), p=length(families$id$binomial))
  }
  if(length(families$id$poisson)>0){
    Y[,families$id$poisson] <- gen_poisson(linpar$linpar[,families$id$poisson,drop=F], n=nrow(linpar$linpar), p=length(families$id$poisson))
  }

  list(Y=Y, Z=Z, linpar=linpar$linpar)
}


# Generate a matrix of covariates
gen_X <- function(nobs, k, intercept){
  X <- matrix(rnorm(nobs*k), nobs, k)
  if(k >= 1 && intercept){
    X[,1] <- 1
  }
  X
}

# Generate a matrix of independent latent factors
gen_Z <- function(nobs, q){
  matrix(stats::rnorm(nobs*q), nobs, q)
}

# Parameters Generation

gen_A <- function(p, q, setting="A", nonzero=100, prop=.4) {
  # setting "A" has only the 100 first loadings that are non-zero
  A <- matrix(runif(p*q, -2, 2), p, q)
  if (setting == "A") {
    if(p>nonzero) {
      A[(nonzero+1):p, ] <- 0
    }
  } else if (setting =="B") {
    shift_size = round(p * (1-prop)/(q-1))

    for (k in 1:q) {

      nonzero_start = (k-1) * shift_size + 1

      if (k == q) {
        nonzero_end = p
      } else {
        nonzero_end   = nonzero_start + round(prop*p) - 1
      }

      nonzeros <- (1:p) %in% (nonzero_start:nonzero_end)
      A[!nonzeros,k] <- 0
    }
  }
  A
}
