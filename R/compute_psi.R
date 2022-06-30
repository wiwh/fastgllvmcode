#' Computes the psi function
#'
#' @param Y: n*p matrix of responses
#' @param Z: n*q matrix of latent variables
#' @param X: n*k design matrix
#' @param A: p*q loadings matrix
#' @param B: p*k coefficients matrix
#' @param family: c("gaussian", "negbin", "poisson", "binomial"), or a named
#' list of indices. See documentation.
#' @param family: known parameters
compute_psi <- function(Y, Z, X, A, B, phi, family, known_par) {

}

# Computing zstar ---------------------

#' Compute zstar
#'
#' Implemented using Newton method. Starting values are important.
#'
#' @inheritParams fastgllvm
#' @param families: a named list of indices corresponding to the families and their respective object
#' @param start: starting values for zstar
#' @param simplify_hessian: if TRUE, only the diagonal of the hessian is used
#'   (not yet implemented)
#' @param line_search: if TRUE, adds a line_search. This makes the convergence
#'   more robust, at the potential cost of computation speed.
#' @param maxit: maximum iterations to compute zstar
#' @param thres: threshold to stop the routine

#TODO: testthat: computes the same first step as glm starting at 0, for all family types.

compute_zstar <- function(Y, A, phi, XB, families, dims, start=NULL, maxit=100, thresh=1e-3, save=F){
  if (is.null(start)) {
    Zstar <- compute_zstar_starting_values(Y, A, XB, families, dims)
  } else {
    stopifnot(is.matrix(start))
    Zstar <- start
  }

  hist <- list()

  # Convergence is independent per i, so we stop computing independently, too
  conv <- rep(F, dims$n)
  crit <- rep(Inf, dims$n)

  linpar <- matrix(NA, dims$n, dims$p)


  Zstar.old <- Zstar
  for(i in 1:maxit){
    # print(beta)
    Zstar.old[!conv,] <- Zstar[!conv, , drop=F]
    if (length(XB) > 0) {
      linpar[!conv,] <- compute_linpar(Zstar[!conv, , drop=F], A, XB=XB[!conv, , drop=F])$linpar
    } else {
      linpar[!conv,] <- compute_linpar(Zstar[!conv, , drop=F], A, XB=XB)$linpar
    }

    Zstar[!conv, ] <- Zstar[!conv, , drop=F] - compute_zstar_Hneg_score(
      compute_zstar_score(Zstar[!conv, , drop=F], Y[!conv, , drop=F], A, phi, linpar[!conv, , drop=F], families),
      compute_zstar_hessian(Y[!conv, , drop=F], A, phi, linpar[!conv, , drop=F], families)
    )

    if (save) hist[[i]] <- as.vector(Zstar)

    crit[!conv] <- rowMeans(abs(Zstar.old[!conv, , drop=F] - Zstar[!conv, , drop=F]))
    conv[!conv] <- crit[!conv] < thresh
    cat("\n", sum(conv), " converged")
    if (all(conv)) break()
  }
  if (save) hist <- do.call(rbind, hist)
  list(Zstar=Zstar, hist=hist, niter=i, converged=ifelse(i==maxit, F, T))
}



compute_zstar_score <- function(Z, Y, A, phi, linpar, families) {
  linpar_bprime <- compute_linpar_bprime(linpar, families)
  (Y - linpar_bprime) %*% (A/phi) - Z
}

compute_zstar_hessian <- function(Y, A, phi, linpar, families) {
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar, families)

  q <- ncol(A)

  if(q > 1) {
    lapply(1:nrow(Y), function(i) {
      - (t(A) %*% (A * linpar_bprimeprime[i,] / phi) + diag(q))
    })
  } else {
    - (t(t(linpar_bprimeprime)*(as.vector(A)/phi)) %*% A + 1)  # TODO: explain this vector trick
  }
}

compute_zstar_Hneg_score <- function(score, hessian){
  if(is.list(hessian)) {
    t(sapply(seq_along(hessian), function(i){
      as.vector(solve(hessian[[i]]) %*% score[i,])
    }))
  } else {
    score/hessian
  }
}

#' Compute Starting values for zstar.
#'
#' Good  starting values are used to speed up / robustify the convergence. They are computed
#' by first transforming the response and then taking the OLS.
compute_zstar_starting_values <- function(Y, A, XB, families, dims) {
  if (length(families$id$poisson) > 0) {
    Y[, families$id$poisson] <- log(Y[, families$id$poisson] + 1e-3)
  }

  if (length(XB) == 0) {
    # no offset
    start <- as.matrix(t(lm(t(Y)~0+A)$coef))
  } else {
    # the offset depends on families...
    offset <- XB
    for (i in seq_along(families$id)) {
      offset[, families$id[[i]]] <- families$object[[i]]$linkinv(XB[,families$id[[i]]])
    }
    start <- as.matrix(t(lm(t(Y)~0+A, offset=t(offset))$coef))

    # This is equivalent to
    # start <- t(sapply(1:dims$n, function(i)as.vector(lm(Y[i,] ~ 0 + A, offset=offset[i,])$coef)))
  }
  colnames(start) <- paste0("Z", 1:dims$q)
  start
}

#' Computes the linear parameter
#'
#' @inheritParams compute_psi
compute_linpar <- function(Z, A, X=NULL, B=NULL, XB=NULL) {
  if(is.null(XB)) {
    if(dim(X)[2] == 0){
      ZA <- Z %*% t(A)
      XB <- vector(length = 0L)
      linpar <- ZA
    } else {
      ZA <- Z %*% t(A)
      XB <- X %*% t(B)
      linpar <- ZA + XB
    }
  } else {
    ZA <- Z %*% t(A)
    if(length(XB) > 0) {
      linpar <- ZA + XB
    } else {
      linpar <- ZA
    }
  }
  list(linpar = linpar, ZA = ZA, XB = XB)
}

compute_linpar_bprime <- function(linpar, families) {
  stopifnot(is.matrix(linpar))
  for(i in seq_along(families$id)){
    linpar[,families$id[[i]]] <- families$objects[[i]]$linkinv(linpar[,families$id[[i]]])
  }
  linpar
}

compute_linpar_bprimeprime <- function(linpar, families) {
  stopifnot(is.matrix(linpar))
  for(i in seq_along(families$id)){
    linpar[,families$id[[i]]] <- families$objects[[i]]$mu.eta(linpar[,families$id[[i]]])
  }
  linpar
}

if(0){
  # TODO: test that this converges
  set.seed(2131)
  f <- gen_fastgllvm(n=1000, p <- 8, q=2, k=4, intercept=T, phi=c(rep(1, p/2), (1:(p/2))/1), family=list(binomial=1:(p/4), poisson=(p/4+1):(p/2), gaussian=(p/2+1):p))
  f <- gen_fastgllvm(n=1000, p <- 10, q=1, k=1, intercept=F, phi=c(rep(1, p/2), (1:(p/2))/1), family=binomial)

  Zstar <- compute_zstar(f$Y, f$A, f$phi, f$linpar$XB, f$families, f$dims, save=F, thresh=1e-3, start=f$Z)
  plot(f$Z, Zstar$Zstar, main = paste("Converged in ", Zstar$niter, " iterations."))

  par(mfrow=c(2,1))
  qqnorm(Zstar$Zstar)
  par(mfrow=c(1,1))
}
