
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

plot_fastgllvm <- function(fastgllvm){
  par(mfrow=c(4,1))
  A.dim <- dim(fastgllvm$hist$A)
  t.seq <- seq(1, A.dim[1], l=min(1000, A.dim[1]))

  sample <- sample(1:A.dim[2], min(250, A.dim[2]))
  ts.plot(fastgllvm$hist$A[t.seq, sample])

  B.dim <- dim(fastgllvm$hist$B)
  sample <- sample(1:B.dim[2], min(250, B.dim[2]))
  ts.plot(fastgllvm$hist$B[t.seq, sample])

  phi.dim <- dim(fastgllvm$hist$phi)
  sample <- sample(1:phi.dim[2], min(250, phi.dim[2]))
  ts.plot(fastgllvm$hist$phi[t.seq, sample])

  ts.plot(fastgllvm$hist$crit[t.seq])
  par(mfrow=c(1,1))
}

#' Print a fastgllvm object
#'
#' @export
print.fastgllvm <- function(fastgllvm){
  cat("test")
}
