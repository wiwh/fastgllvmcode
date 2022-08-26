
# Methods
# -------

#' Predict the latent variable Z from a gllvm, assuming it is multivariate normally distributed.
#'
#' @param f: an object of class "fastgllvm"
#' @param method: one of "glm" or "ridge" (default)
#'
#' @return a n x q matrix of factors
#'
#' @export
predict.fastgllvm <- function(f, method=c("fastgllvm", "glm", "glmnet")){
  compute_zstar(f$Y, f$A, f$linpar$XB, f$families)
}

# Plot the fastgllvm object
#' @export
plot.fastgllvm <- function(f){
  if(length(f$fit) == 0) {
    stop("Fit the model before attempting to plot.")
  }
  par(mfrow=c(3,1))
  ts.plot(f$fit$hist$A[,1:min(100, f$dimensions$p*f$dimensions$q)], main="Convergence plot: loadings.", xlab="Iteration", col=1:min(100, f$dimensions$p*f$dimensions$q), lwd=2)
  if(!is.null(f$fit$hist$B)) {
    ts.plot(f$fit$hist$B[,1:min(100, f$dimensions$p)], main="Convergence plot: betas.", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  }
  ts.plot(f$fit$hist$phi[,1:min(100, f$dimensions$p)], main = "Convergence plot: communalities.", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  par(mfrow=c(1,1))
}

#' Update the fit with another round of stochastic approximation
#' @export
update.fastgllvm <- function(f, ...){
  if(length(f$fit) == 0) {
    stop("Fit the model before attempting to update it.")
  }

  arguments <- list(...)

  if (is.null(arguments[["controls"]])) {
    controls <- f$fit$controls
    controls$alpha <- controls$alpha/2
  } else {
    controls <- arguments[["controls"]]
  }

  if (is.null(arguments[["verbose"]])) {
    verbose <- F
  }

  if (is.null(arguments[["hist"]])) {
    hist <- T
  }
  fastgllvm.fit(f, parameters.init = f$parameters, controls=controls, verbose=verbose, hist=hist)
}

simulate_fastgllvm <- function(fastgllvm, n=NULL){
  nsim <- ifelse(is.null(n), fit$n, n)
  with(fastgllvm, {
    Z <- gen_Z(n, q)
    if(nsim!=n){
      sampl <- sample(n, nsim, replace=T)
    } else {
      sampl <-1:n
    }
    gen_Y(A, B, phi, Z[sampl,, drop=F], X[sampl,, drop=F], family)
  })
}


#' Print a fastgllvm object
#'
#' @export
print.fastgllvm <- function(fastgllvm){
  cat("The 'print' method for a 'fastgllvm' object has not been implemented yet.")
}
