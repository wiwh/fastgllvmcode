
# Methods
# -------

#' Predict the latent variable Z from a gllvm, assuming it is multivariate normally distributed.
#'
#' @param f: an object of class "gllvmprime"
#' @param method: one of "glm" or "ridge" (default)
#'
#' @return a n x q matrix of factors
#'
#' @export
predict.gllvmprime <- function(f, method=c("gllvmprime", "glm", "glmnet")){
  compute_zstar(f$Y, f$A, f$linpar$XB, f$families)
}

# Plot the gllvmprime object
#' @export
plot.gllvmprime <- function(f, plot.last=NULL, plot.covZ=F){
  if(is.null(f$hist)) {
    stop("Fit the model before attempting to plot.")
  }
  plot_id <- 1:min(100, f$dimensions$p*f$dimensions$q)
  iter <- (1:nrow(f$hist$A))
  if(!is.null(plot.last))  iter <- iter[(length(iter) - min(length(iter), plot.last) + 1):length(iter)]

  if (plot.covZ) {
    par(mfrow=c(3,2))
  } else {
    par(mfrow=c(2,2))
  }
  ts.plot(f$hist$A[iter,plot_id], main="Convergence plot: loadings", xlab="Iteration", col=plot_id, lwd=2)
  points(rep(iter[length(iter)], length(plot_id)), as.vector(f$parameters$A)[plot_id], col=plot_id)
  if(!is.null(f$hist$B)) {
    ts.plot(f$hist$B[iter,1:min(100, f$dimensions$p)], main="Convergence plot: betas", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  }
  ts.plot(f$hist$phi[iter,1:min(100, f$dimensions$p)], main = "Convergence plot: scales", xlab="Iteration", col=1:min(100, f$dimensions$p), lwd=2)
  ts.plot(f$hist$deviance[iter], main = "Convergence plot: deviance", xlab="Iteration", lwd=2, ylab="")
  if (plot.covZ) {
    ts.plot(f$hist$covZ, main = "Convergence plot: covZ", xlab="Iteration", col=1:f$dimensions$q**2, lwd=2)
  }
  par(mfrow=c(1,1))
}

#' Update the fit with another round of stochastic approximation
#' @param ...: named arguments, same as in the call to gllvmprime
#' @method update gllvmprime
#' @export
"update.gllvmprime" <- function(fg, ...){
  args <- list(...)

  # unpacking control arguments
  for (name in names(fg$controls)) {
    if (!is.null(args[[name]])) {
      fg$controls[[name]] <- args[[name]]
    }
  }

  gllvmprime.fit(fg, parameters.init = fg$parameters, controls=fg$controls)
}

#' Subsets the gllvmprime object to only include observational units indicated in `index`
#'
#' @param index: numeric vector indicating which observations to keep
#' @export
subset.gllvmprime  <- function(gllvmprime, index) {
  if (is.null(index)) {
    index <- 1:gllvmprime$dimensions$n
  }
  gllvmprime$dimensions$n <- length(index)
  gllvmprime$Z <- gllvmprime$Z[index,,drop=F]
  gllvmprime$Y <- gllvmprime$Y[index,,drop=F]
  gllvmprime$X <- gllvmprime$X[index,,drop=F]
  if (!is.null(gllvmprime$Miss))gllvmprime$Miss <- gllvmprime$Miss[index,,drop=F]
  if (!is.null(gllvmprime$linpar)) gllvmprime$linpar <- gllvmprime$linpar[index,,drop=F]
  if (!is.null(gllvmprime$mean)) gllvmprime$mean <- gllvmprime$mean[index,,drop=F]
  gllvmprime
}


#' Simulate new data from a given gllvmprime object
#´
#' @param index: numeric vector indicating which observation unit to simulate, useful when covariate information X is provided.
#' @param value: a gllvmprime object with the simulations
#' @param conditional: should the simulation be conditional on gllvmprime$Z?
#' @export
simulate.gllvmprime <- function(gllvmprime, nsim=1, conditional=F, return_object=F){
  if(!conditional) {
    gllvmprime$Z <- with(gllvmprime, gen_Z(dimensions$n, dimensions$q))
  }
  simu <- with(gllvmprime, gen_Y(Z=Z, X=X, parameters = parameters, families = families))
  if (return_object) {
    gllvmprime$Y <- simu$Y
    gllvmprime$Z <- simu$Z
    gllvmprime$linpar <- simu$linpar
    gllvmprime$mean <- NULL
    gllvmprime$Miss <- NULL
    return(gllvmprime)
  } else {
    return(simu)
  }
}


#' Print a gllvmprime object
#'
#' @export
print.gllvmprime <- function(gp, n=10){
  if (is.null(gp$fit)) {
    cat("\nThe model has not been fit to data yet. Run `gllvmprime` on data.")
  } else {
    cat("Fitted GLLVM model:")
    cat("\n-------------------")
    A <- gp$parameters$A
    colnames(A) <- NULL
    B <- gp$parameters$B
    colnames(B) <- NULL
    A <- dplyr::tibble(A=A)
    B <- dplyr::tibble(B=B)
    phi <- gp$parameters$phi
    cat("\nLoadings:\n")
    print(A, n=n)
    cat("\nFixed coefficients:\n")
    print(B, n=n)
    cat("\nScale parameters:\n")
    print(phi)
    cat("\nNumerical controls:")
      cat("\n\tNumber of iterations: \tmaxit =", gp$controls$maxit)
      cat("\n\tBatch size: \t\tbatch_size =", gp$controls$batch_size)
      cat("\n\tEnd value for alpha: \talpha =", gp$controls$alpha)
  }
}
