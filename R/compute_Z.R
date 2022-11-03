#' Compute the MAP
#'
#' @description For each observational unit i, computes using the Newton method the maximum a posteriori value, i.e. the values that maximizes the distribution of Z|Y_i for given parameter values.
#'
#' @inheritParams fastgllvm
#' @param start: initial estimate of the MAP. If NULL, initial values will be computed depending based on an initial estimate of mu.
#' @param maxit: maximum number of Newton iterations
#' @param thresh: threshold to stop the routine; see
#' @param lambda: the regularization parameter
#' @param save: iterations are saved
#' @param verbose: if TRUE, documents and prints the progress
#'
#' @details
#' For computational speed, the convergence is checked for each observational unit `i` separately, and computations are vectorized across observations as much as possible.
#' The `lambda` parameter is the regularization parameter. Setting it to 1 (the fault) returns the MAP. Higher values regularizes the values towards 0.

compute_Z <- function(fg, start=NULL, maxit=10, thresh=1e-3, lambda=1, save=F, verbose=F){
  with(fg,{
    with(parameters, {
      if(is.null(B)) {
        XB <- NULL
      } else {
        XB <- X %*% t(B)
      }

      if (is.null(start)) {
        Z <- compute_Z_starting_values(Y, A, XB, families)
      } else {
        stopifnot(is.matrix(start))
        Z <- start
      }

      hist <- list()

      # Initialize convergence vectors for all rows
      conv <- rep(F, nrow(Y))
      crit <- rep(Inf, nrow(Y))

      # Initialize the linear parameter
      linpar <- linpar_bprime <- matrix(NA, nrow(Y), ncol(Y))


      Z.old <- Z

      for(i in 1:maxit){
        Z.old[!conv,] <- Z[!conv, , drop=F]
        if (length(XB) > 0) {
          linpar[!conv,] <- compute_linpar(Z[!conv, , drop=F], A, XB=XB[!conv, , drop=F])$linpar
        } else {
          linpar[!conv,] <- compute_linpar(Z[!conv, , drop=F], A, XB=XB)$linpar
        }

        linpar_bprime[!conv,] <- compute_linpar_bprime(linpar[!conv,,drop=F], families)

        Z[!conv, ] <- Z[!conv, , drop=F] - Z_update(
          compute_Z_score(Z[!conv, , drop=F], Y[!conv, , drop=F], A, phi, linpar_bprime[!conv, , drop=F], families, lambda=lambda),
          compute_Z_hessian(A, phi, linpar[!conv, , drop=F], families, lambda=lambda)
        )

        if (save) hist[[i]] <- as.vector(Z)

        crit[!conv] <- rowMeans(abs(Z.old[!conv, , drop=F] - Z[!conv, , drop=F]))
        conv[!conv] <- crit[!conv] < thresh
        if (verbose) cat("\n", sum(conv), " converged")
        if (all(conv)) break()
      }
      if (save) hist <- do.call(rbind, hist)
      list(Z=Z, linpar=linpar, linpar_bprime=linpar_bprime, hist=hist, niter=i, converged=ifelse(i==maxit, F, T))
    })
  })
}




compute_Z_score <- function(Z, Y, A, phi, linpar_bprime, families, lambda) {
  # For justification of this computation, check the Rmd file "NA".
  # We simply ignore the missing values in the score: careful the Hessian
  # needs to have the same simplification, lest the scores be divided by the wrong quantity:
  # there need to be the same elements ignored in the scores as they are in the Hessian.

  # if (!is.null(Miss)) {
  #   Y[Miss] <- 0 #TODO: do this in parent!
  #   linpar_bprime[Miss] <- 0
  # }

  (Y - linpar_bprime) %*% (A/phi) - Z*lambda
}

#TODO: we need A/phi here as well: it will never change for Z.
compute_Z_hessian <- function(A, phi, linpar, families, lambda) {
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar, families)
  q <- ncol(A)

  # TODO: this should be done with tensor products
  if(q > 1) {
    lapply(1:nrow(linpar_bprimeprime), function(j) {
      - (t(A) %*% (A * (linpar_bprimeprime[j,] / phi)) + diag(q)*lambda)
    })
  } else {
    - (t(t(linpar_bprimeprime)*(as.vector(A)/phi)) %*% A + lambda)  # TODO: explain this vector trick
  }
}

# Compute the update for Z, as the negative of the hessian times the score.
# If the hessian is a vector (method == "newton"), the computations are made row-wise. Else they are done as an element-wise matrix multiplication (method=="quasi").
Z_update <- function(score, hessian){
  if(is.list(hessian)) {
    t(sapply(seq_along(hessian), function(i){
      solve(hessian[[i]], score[i,])
    }))
  } else {
    score/hessian
  }
}


#' Compute Starting values for Z.
#'
#' Good  starting values are used to speed up / robustify the convergence. They are computed
#' by first transforming the response and then taking the OLS.
compute_Z_starting_values <- function(Y, A, XB, families) {
  if (length(families$id$poisson) > 0) {
    Y[, families$id$poisson] <- log(Y[, families$id$poisson] + 1e-3)
  }

  if (is.null(XB)) {
    # no offset
    start <- as.matrix(t(lm(t(Y)~0+A)$coef))
  } else {
    # TODO: making offset dependent of families (as we would think) DOES NOT WORK WELL for initialization when k is big.
    # UNDERSTAND WHY
    # DO NOT CHANGE OFFSET, but if you want, this is how: for (i in seq_along(families$id)) {
    # DO NOT CHANGE OFFSET, but if you want, this is how:   offset[, families$id[[i]]] <- families$object[[i]]$linkinv(XB[,families$id[[i]]])
    # DO NOT CHANGE OFFSET, but if you want, this is how: }

    start <- as.matrix(t(lm(t(Y)~0+A, offset=t(XB))$coef))
    # This is equivalent, but faster, to
    # start <- t(sapply(1:nrow(Y), function(i)as.vector(lm(Y[i,] ~ 0 + A, offset=offset[i,])$coef)))
    # Alternatively, it is also equivalent to
    # start <- (Y - offset) %*% (A %*% solve(t(A) %*% A))
  }

  colnames(start) <- paste0("Z", 1:ncol(A))
  start
}


if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 100
  gaussian <- 100
  binomial <- 100
  nobs <- 100
  q <- 2
  p <- poisson + gaussian + binomial

  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(10030)
  fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, phi=runif(p) + 0.5, miss.prob = 0, scale=1)
  zhat <- compute_Z(fg, start=fg$parameters$Z, lambda=1)

  plot(fg$Z, zhat$Z, xlim=c(-3,3), ylim=c(-3,3)); abline(0, 1, col=2)
}
