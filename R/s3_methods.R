
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

fit_fastgllvm <- function(fastgllvm, method="SA", H=1, maxit=250, tol=1e-5, learning_rate = NULL,  learning_rate.args = NULL, verbose = T ){
  if(is.null(learning_rate)) learning_rate <- ifelse(method=="SA", "exp", "constant")
  if(is.character(learning_rate)){
    learning_rate <- ff_learning_rate(method=learning_rate, maxit=maxit, learning_rate.args = learning_rate.args)
  }
  learning_rate.seq <- learning_rate(1:maxit)
  cat(nrow(fastgllvm$hist$A))
  class(fastgllvm) <- NULL
  fastgllvm <- within(fastgllvm,
  {
    hist.i <- nrow(hist$A)
    cat(" - ", hist.i, " \n")

    hist$A <- rbind(hist$A, matrix(0, maxit, p*q))
    hist$B <- rbind(hist$B, matrix(0, maxit, p*k))
    hist$phi <- rbind(hist$phi, matrix(0, maxit, p))
    hist$crit <- c(hist$crit, rep(0, maxit))

    ################################################################################################################
    Y.c <- scale(Y, scale=F)  # TODO CHANGE THIS BACK TODO
    generate_Z <- generate_Z_functionfactory(n, q, method=method, H=H)

    i <- 0
    converged <- FALSE
    while(i < maxit & !converged){
      i <- i+1
      Psi <- get_Psi(Y, Y.c, A, B, phi, X, family, generate_Z)

      # udate A
      A <- A + learning_rate.seq[i] * Psi$A
      broken.A <- abs(A)>10
      A[broken.A] <- 10 * sign(A[broken.A])

      # update B
      B  <- B + learning_rate.seq[i] * Psi$B

      # update phi
      phi <- phi + learning_rate.seq[i] * Psi$phi

      # save
      hist$A[hist.i + i, ] <- as.vector(A)
      hist$B[hist.i + i, ] <- as.vector(B)
      hist$phi[hist.i + i, ] <- as.vector(phi)

      hist$crit[hist.i + i] <- learning_rate.seq[i] * norm(Psi$A)/(p*q)
      if(hist$crit[hist.i + i] < tol) converged <- TRUE

      if(verbose)cat("\ni: ", i, " - norm:", hist$crit[hist.i + i], " - learning rate:", learning_rate.seq[i])
      # check if the criterion is small enough to jump to the next "repetition", where the learning rate increases again
      if(converged){
        # fill in the histories
        hist$A <- hist$A[1:(hist.i + i),]
        hist$B <- hist$B[1:(hist.i + i),]
        hist$phi <- hist$phi[1:(hist.i + i),]
        hist$crit <- hist$crit[1:(hist.i + i)]
      }
    }
  })
  class(fastgllvm) <- "fastgllvm"
  fastgllvm
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
