# Predict z as the MAP.
# returns the fitted means, too
bernoulli.predict <- function(Y, par, X=NULL, maxiter=1000, eps=1e-5, verbose=F){
  n <- nrow(Y)

  if(!is.null(X)){
    fixedPart <- X %*% par$B
  } else {
    fixedPart <- 0
  }

  # initialize Latent
  Z <- matrix(0, n, par$q)

  conv <- FALSE
  iter <- 1
  while(!conv & iter < maxiter){
    Z.old <- Z
    iter <- iter + 1

    latPart <- Z %*% t(par$A)
    linpar <- latPart + fixedPart
    Z <- Z + .05*(-Z + (Y - sigmoid(linpar))%*%par$A)
    crit <- mean((Z-Z.old)^2)
    conv <- crit < eps
    if(verbose)cat("\niter:", iter, crit, " conv ", conv)
  }

  list(linpar=linpar, Z=Z, converged=conv)
}
