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
compute_psi <- function (Y, Z, X, A, B, phi, families) {
}

compute_pi <- function(Y, Z, X, A, B, phi, families, maxit=100) {
  # TODO: this can easily be palalellized
  sol <- t(sapply(1:ncol(Y), function(j){
    compute_glm(j, Y, Z, X, A, B, phi, families, maxit=maxit)
  }))
  A <- sol[, startsWith(colnames(sol), "A"), drop=F]
  B <- sol[, startsWith(colnames(sol), "B"), drop=F]
  phi <- sol[, startsWith(colnames(sol), "phi")]

  list(A=A, B=B, phi=phi)
}

compute_glm <- function(j, Y, Z, X, A, B, phi, families, maxit=100) {
  if(is.null(A) | is.null(B) | is.null(phi)) {
    fit <- glm(Y[,j]~0+X+Z, family=families$vec[j], maxit=maxit)
  } else {
    fit <- glm(Y[,j]~0+X+Z, start=c(B[j,], A[j,]), family=families$vec[j], maxit=maxit)
  }
  B_j <- fit$coef[1:ncol(X)]
  A_j <- fit$coef[(ncol(X)+1):length(fit$coef)]
  phi_j <- 1
  if (family=="gaussian") {
    phi <- var(fit$residuals)
  }
  c(A_j=A_j, B_j=B_j, phi_j=phi_j)
}

compute_psi_AB <- function (Y, Z, X, B, A, phi, families, lambda=0) {
  linpar <- compute_linpar(Z, A, X=X, B=B)
  linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)

  psi_A <- do.call(rbind, lapply(1:ncol(Y), function(j){
    as.vector(solve(
      compute_psi_A_j_hessian(Z, phi[j], linpar_bprimeprime[,j]) - diag(lambda, ncol(Z)),
      compute_psi_A_j(Y[,j], Z, phi[j], linpar_bprime[,j]) - lambda * A[j,]
    ))
  }))

  # psi_B <- do.call(rbind, lapply(1:ncol(Y), function(j){
  #   as.vector(solve(
  #     compute_psi_A_j_hessian(X, phi[j], linpar_bprimeprime[,j]) + diag(lambda, ncol(Z)),
  #     compute_psi_A_j(Y[,j], X, phi[j], linpar_bprime[,j]) + lambda * B[j,]
  #   ))
  # }))

  # psi_AB <- do.call(rbind, lapply(1:ncol(Y), function(j){
  #   as.vector(solve(
  #     compute_psi_A_j_hessian(cbind(X, Z), phi[j], linpar_bprimeprime[,j]),
  #     compute_psi_A_j(Y[,j], cbind(X,Z), phi[j], linpar_bprime[,j])
  #   ))
  # }))

  list(psi_A=psi_A) #, psi_B=psi_B) #, psi_AB=psi_AB)
}

compute_psi_A_j <- function(Y_j, Z, phi_j, linpar_bprime_j){
  (t(Z) %*% (Y_j - linpar_bprime_j))/phi_j
}

compute_psi_A_j_hessian <- function(Z, phi_j, linpar_bprimeprime_j){
  -(t(Z) %*% (Z*(linpar_bprimeprime_j)))/phi_j
}

compute_psi_B <- function () {

}


compute_psi_phi <- function () {

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

compute_zstar <- function(Y, A, phi, XB, families, dims, start=NULL, maxit=100, thresh=1e-3, save=F, scale=F){
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

    if (scale) Zstar <- scale_zstar(Zstar)
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

scale_zstar <- function(Z){
  nq <- prod(dim(Z))
  gaussian_quantiles <- qnorm((1:nq)/(nq+1)) # TODO: don't run that each time...
  Z[order(Z)] <- gaussian_quantiles
  Z
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

