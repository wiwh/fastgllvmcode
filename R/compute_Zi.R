
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
#' @param Y_i: a vector of observations
#' @param A: the loadings (act as the covariates)
#' @param XB_i: the offset
#' @param Miss_i: the missing values

#TODO: testthat: computes the same first step as glm starting at 0, for all family types.

compute_zstar_i <- function(Y_i, A, phi, XB_i, families, start=NULL, maxit=100, thresh=1e-3, save=F, scale=F, Miss_i=NULL){
  if (is.null(start)) {
    Zstar <- compute_zstar_starting_values(Y_i, A, XB_i, families, Miss_i)
  } else {
    stopifnot(is.matrix(start))
    Zstar <- start
  }

  hist <- list()
  linpar <- matrix(NA, ncol(Y))


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



compute_zstar_i_score <- function(Z_i, Y_i, A, phi, linpar_i, families, Miss_i) {
  linpar_bprime_i <- compute_linpar_bprime(linpar_i, families)

  if (is.null(Miss)) {

    score <- as.vector(t(A/phi) %*% (Y_i - linpar_bprime_i))  - Z_i

  } else {

    score <- as.vector(t(A[!Miss_i,]/phi[!Miss_i]) %*% (Y_i[!Miss_i] - linpar_bprime_i[!Miss_i]))  - Z_i

  }
  score
}

compute_zstar_i_hessian <- function(Y, A, phi, linpar, families) {
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

compute_zstar_i_Hneg_score <- function(score, hessian){
  if(is.list(hessian)) {
    t(sapply(seq_along(hessian), function(i){
      as.vector(solve(hessian[[i]]) %*% score[i,])
    }))
  } else {
    score/hessian
  }
}

scale_zstar_i <- function(Z){
  nq <- prod(dim(Z))
  gaussian_quantiles <- qnorm((1:nq)/(nq+1)) # TODO: don't run that each time...
  Z[order(Z)] <- gaussian_quantiles
  Z
}

#' Compute Starting values for zstar.
#'
#' Good  starting values are used to speed up / robustify the convergence. They are computed
#' by first transforming the response and then taking the OLS.
#'
#' NAs are automatically taken care of by lm.

compute_zstar_i_starting_values <- function(Y_i, A, XB_i, families) {
  if (length(families$id$poisson) > 0) {
    Y_i[families$id$poisson] <- log(Y_i[families$id$poisson] + 1e-3)
  }

  if (is.null(XB_i)) {
    # no offset
    start <- as.vector(lm(Y_i~0+A, na.action = na.omit)$coef)
  } else {
    # the offset depends on families...
    offset <- XB_i
    for (i in seq_along(families$id)) {
      offset[families$id[[i]]] <- families$object[[i]]$linkinv(XB_i[families$id[[i]]])
    }
    offset[is.na(Y_i)] <- NA
    start <- as.vector(lm(Y_i~0+A, offset=offset, na.action = na.omit)$coef)
  }

  names(start) <- paste0("Z", 1:ncol(A))
  start
}




if(0) {
  fg <- gen_fastgllvm(nobs=1000, p=40, q=1, family="binomial", k=0, intercept=F)
  zhat <- compute_zstar(fg$Y, fg$parameters$A, fg$parameters$phi, fg$linpar$XB, fg$families, Miss=NULL)
  zhat_i <- compute_zstar(fg$Y[1,], fg$parameters$A, fg$parameters$phi, fg$linpar$XB[1,], fg$families)
  plot(fg$Z, zhat$Zstar)

  fg_miss <- gen_fastgllvm(nobs=1000, p=40, q=1, family="binomial", k=1, intercept=T, miss.prob = 0.1)
  zhat <- compute_zstar_starting_values(fg_miss$Y, fg_miss$parameters$A, fg_miss$linpar$XB, fg_miss$families, fg_miss$Miss)
  plot(fg_miss$Z, zhat)


  library(mirtjml)
  res <- mirtjml::mirtjml_expr(fg$Y, fg$dimensions$q, cc=1000, tol = 1e-2)
  res <- mirtjml::mirtjml_expr(fg$Y, fg$dimensions$q)
  plot(res$A_hat, fg$parameters$A); abline(0,1,col=2)
  zhat <- compute_zstar(fg$Y, res$A_hat, rep(1, fg$dimensions$p), res)
  plot(fg$parameters$A, res$A_hat)
  plot(fg$Z, res$theta_hat)
  abline(0,1,col=2); abline(0,-1,col=2)
}
