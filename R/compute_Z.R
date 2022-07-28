
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

compute_zstar <- function(Y, A, phi, XB, families, start=NULL, maxit=100, thresh=1e-3, save=F, Miss=NULL, verbose=F){
  if (is.null(start)) {
    Zstar <- compute_zstar_starting_values(Y, A, XB, families, Miss)
  } else {
    stopifnot(is.matrix(start))
    Zstar <- start
  }

  hist <- list()

  # The convergence of our algorithm is row-wise. If a row has converged, we remove it from updating in the next iteration.
  conv <- rep(F, nrow(Y))
  crit <- rep(Inf, nrow(Y))

  linpar <- matrix(NA, nrow(Y), ncol(Y))


  Zstar.old <- Zstar

  # TODO: This next part can be split into groups for parallel computing.
  for(i in 1:maxit){
    Zstar.old[!conv,] <- Zstar[!conv, , drop=F]
    if (length(XB) > 0) {
      linpar[!conv,] <- compute_linpar(Zstar[!conv, , drop=F], A, XB=XB[!conv, , drop=F])$linpar
    } else {
      linpar[!conv,] <- compute_linpar(Zstar[!conv, , drop=F], A, XB=XB)$linpar
    }

    Zstar[!conv, ] <- Zstar[!conv, , drop=F] - compute_zstar_Hneg_score(
      compute_zstar_score(Zstar[!conv, , drop=F], Y[!conv, , drop=F], A, phi, linpar[!conv, , drop=F], families, Miss[!conv, , drop=F]),
      compute_zstar_hessian(A, phi, linpar[!conv, , drop=F], families, Miss[!conv, , drop=F])
    )

    if (save) hist[[i]] <- as.vector(Zstar)

    crit[!conv] <- rowMeans(abs(Zstar.old[!conv, , drop=F] - Zstar[!conv, , drop=F]))
    conv[!conv] <- crit[!conv] < thresh
    if(verbose) cat("\n", sum(conv), " converged")
    if (all(conv)) break()
  }
  if (save) hist <- do.call(rbind, hist)
  list(Zstar=Zstar, hist=hist, niter=i, converged=ifelse(i==maxit, F, T))
}



compute_zstar_score <- function(Z, Y, A, phi, linpar, families, Miss) {
  linpar_bprime <- compute_linpar_bprime(linpar, families)

  # For justification of this computation, check the Rmd file "NA".
  # We simply ignore the missing values in the score: careful the Hessian
  # needs to have the same simplification, lest the scores be divided by the wrong quantity:
  # there need to be the same elements ignored in the scores as they are in the Hessian.

  if (!is.null(Miss)) {
    Y[Miss] <- 0 #TODO: do this in parent!
    linpar_bprime[Miss] <- 0
  }

  (Y - linpar_bprime) %*% (A/phi) - Z
}

#TODO: we need A/phi here as well: it will never change for zstar.
compute_zstar_hessian <- function(A, phi, linpar, families, Miss) {
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar, families)

  q <- ncol(A)

  if(!is.null(Miss)) {
    # TODO: do this in parent!
    linpar_bprimeprime[Miss] <- 0
  }

  # TODO: this should be done with tensor products
  if(q > 1) {
    lapply(1:nrow(linpar_bprimeprime), function(j) {
      - (t(A) %*% (A * linpar_bprimeprime[j,] / phi) + diag(q))
    })
  } else {
    - (t(t(linpar_bprimeprime)*(as.vector(A)/phi)) %*% A + 1)  # TODO: explain this vector trick
  }
}

compute_zstar_Hneg_score <- function(score, hessian){
  if(is.list(hessian)) {
    t(sapply(seq_along(hessian), function(i){
      solve(hessian[[i]], score[i,]) #TODO: regularization should be added as an option
    }))
  } else {
    score/hessian
  }
}


# not used anymore (should not be used here, but maybe in parent... A should be rescaled appropriately)
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
compute_zstar_starting_values <- function(Y, A, XB, families, Miss) {

  if (length(families$id$poisson) > 0) {
    Y[, families$id$poisson] <- log(Y[, families$id$poisson] + 1e-3)
  }

  if (is.null(XB)) {
    # no offset
    if (is.null(Miss)) {
      start <- as.matrix(t(lm(t(Y)~0+A)$coef))
    } else {
      start <- t(sapply(1:nrow(Y), function(i)as.vector(lm(Y[i,] ~ 0 + A, na.action = na.omit)$coef)))
    }
  } else {
    # the offset depends on families...
    offset <- XB
    # TODO: making offset dependent of families (as we would think) DOES NOT WORK WELL for initialization when k is big.
    # UNDERSTAND WHY
    # DO NOT CHANGE OFFSET, but if you want, this is how: for (i in seq_along(families$id)) {
    # DO NOT CHANGE OFFSET, but if you want, this is how:   offset[, families$id[[i]]] <- families$object[[i]]$linkinv(XB[,families$id[[i]]])
    # DO NOT CHANGE OFFSET, but if you want, this is how: }

    if (is.null(Miss)) {
      start <- as.matrix(t(lm(t(Y)~0+A, offset=t(offset))$coef))
      # This is equivalent, but faster, to
      # start <- t(sapply(1:nrow(Y), function(i)as.vector(lm(Y[i,] ~ 0 + A, offset=offset[i,])$coef)))
      # Alternatively, it is also equivalent to
      # start <- (Y - offset) %*% (A %*% solve(t(A) %*% A))

    } else {
    # If there are missing values, compute them individually, lm takes care to omit them, but needs to do so individually:
      start <- t(sapply(1:nrow(Y), function(i)as.vector(lm(Y[i,] ~ 0 + A, offset=offset[i,], na.action = na.omit)$coef)))
    }
  }

  colnames(start) <- paste0("Z", 1:ncol(A))
  start
}

if(0) {
  devtools::load_all()
  set.seed(1231)
  fg <- gen_fastgllvm(nobs=1000, p=1000, q=5, family=c(rep("poisson", 950), rep("gaussian", 20), rep("binomial",30)), k=0, intercept=F, miss.prob = 0.1)

  parameters.init <- initialize_parameters(fg, target=fg$parameters$A)
  if(!is.null(parameters.init$B)) {
    parameters.init$XB <- fg$X %*% t(parameters.init$B)
  } else {
    parameters.init$XB <- NULL
  }
  zhat <- compute_zstar(fg$Y, parameters.init$A, parameters.init$phi, parameters.init$XB, fg$families, Miss=fg$Miss)
  plot(fg$parameters$A, parameters.init$A)
  points(fg$Z, zhat$Zstar, col=2); abline(0,1,col=3)
  # now we rescale.

  # zstart must be the same individual or vector
  zstart <- with(fg, compute_zstar_starting_values(Y, parameters$A, linpar$XB, families, Miss=Miss))

  zstart_individual <- with(fg, {
      t(sapply(1:nrow(Y), function(i){
        compute_zstar_i_starting_values(Y[i,], parameters$A, linpar$XB[i,], families)
      }))
  })

  all.equal(zstart, zstart_individual)

  zstart <- with(fg_miss, compute_zstar_starting_values(Y, parameters$A, linpar$XB, families))
  zstart2 <- with(fg_miss, {
      t(sapply(1:nrow(Y), function(i){
        compute_zstar_i_starting_values(Y[i,], parameters$A, linpar$XB[i,], families)
      }))
  })

  all.equal(zstart, zstart2)
}
