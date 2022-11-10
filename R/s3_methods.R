
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
  plot_id <- 1:min(100, f$dimensions$p*f$dimensions$q)
  par(mfrow=c(2,2))
  ts.plot(f$fit$hist$A[,plot_id], main="Convergence plot: loadings.", xlab="Iteration", col=plot_id, lwd=2)
  points(rep(nrow(f$fit$hist$A), length(plot_id)), as.vector(f$parameters$A)[plot_id], col=plot_id)
  if(!is.null(f$fit$hist$B)) {
    ts.plot(f$fit$hist$B[,1:min(100, f$dimensions$p)], main="Convergence plot: betas.", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  }
  ts.plot(f$fit$hist$phi[,1:min(100, f$dimensions$p)], main = "Convergence plot: communalities.", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  ts.plot(f$fit$hist$covZ, main = "Convergence plot: covZ.", xlab="Iteration", col=1:f$dimensions$q**2, lwd=2)
  par(mfrow=c(1,1))
}

#' Update the fit with another round of stochastic approximation
#' @param ...: named arguments, same as in the call to fastgllvm
#' @method update fastgllvm
#' @export
"update.fastgllvm" <- function(fg, ...){
  args <- list(...)

  # unpacking the extra args
  for (name in names(fg$controls)) {
    if (!is.null(args[[name]])) {
      fg$controls[[name]] <- args[[name]]
    }
  }

  fastgllvm.fit(fg, parameters.init = fg$parameters, controls=fg$controls)
}

#' Subsets the fastgllvm object to only include observational units indicated in `index`
#'
#' @param index: numeric vector indicating which observations to keep
#' @export
subset.fastgllvm  <- function(fastgllvm, index) {
  if (is.null(index)) {
    index <- 1:fastgllvm$dimensions$n
  }
  fastgllvm$dimensions$n <- length(index)
  fastgllvm$Z <- fastgllvm$Z[index,,drop=F]
  fastgllvm$Y <- fastgllvm$Y[index,,drop=F]
  fastgllvm$X <- fastgllvm$X[index,,drop=F]
  if (!is.null(fastgllvm$linpar)) fastgllvm$linpar <- fastgllvm$linpar[index,,drop=F]
  if (!is.null(fastgllvm$mean)) fastgllvm$mean <- fastgllvm$mean[index,,drop=F]
  fastgllvm
}


#' Simulate new data from a given fastgllvm object
#´
#' @param index: numeric vector indicating which observation unit to simulate, useful when covariate information X is provided.
#' @param value: a fastgllvm object with the simulations
#' @param conditional: should the simulation be conditional on fastgllvm$Z?
#' @export
simulate.fastgllvm <- function(fastgllvm, nsim=1, conditional=F, return_fastgllvm=F){
  if(!conditional) {
    fastgllvm$Z <- with(fastgllvm, gen_Z(dimensions$n, dimensions$q))
  }
  simu <- with(fastgllvm, gen_Y(Z=Z, X=X, parameters = parameters, families = families))
  if (return_fastgllvm) {
    fastgllvm$Y <- simu$Y
    fastgllvm$Z <- simu$Z
    fastgllvm$linpar <- simu$linpar
    fastgllvm$mean <- NULL
    fastgllvm$Miss <- NULL
    return(fastgllvm)
  } else {
    return(simu)
  }
}


#' Print a fastgllvm object
#'
#' @export
print.fastgllvm <- function(fastgllvm){
  cat("The 'print' method for a 'fastgllvm' object has not been implemented yet.")
}
