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



predict_z <- function(Y, A, B, phi, X, family){
  offset <- X %*% t(par$B)
  z <- t(sapply(1:nrow(Y), function(i) glm(Y[i,]~ 0 + par$A, family=family, offset = offset[i,])$coef))
  z
}

if(0){
  set.seed(2145)
  family <- binomial()
  g <- gen_fastgllvm(n = 100, p = 1000, q = 1, family =family, k=5)
  # g$par$B <- g$par$B * 0
  z <- with(g, predict_z(Y, A, B, phi, X, family))
  plot(z, g$Z); abline(0,1,col=2)
}
