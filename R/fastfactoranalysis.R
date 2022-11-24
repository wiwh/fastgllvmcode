# gen_data

ffa_gen_Z <- function(n, q) {
  matrix(rnorm(n*q), n, q)
}

ffa_gen_Y <- function(n=NULL, p=p, q=q, A=NULL, Psi=NULL, Z=NULL) {
  if(is.null(A)) A <- ffa_gen_L(p, q)
  p <- nrow(A)
  q <- ncol(A)
  if(is.null(Z)) Z <- ffa_gen_Z(n, q)
  n <- nrow(Z)
  if(is.null(Psi)) Psi <- ffa_gen_Psi(p)

  linpar <- Z %*% t(A)

  Y <- sapply(1:p, function(j) {
    rnorm(n, linpar[,j], sd=sqrt(Psi[j]))
  })
  list(Y=Y, A=A, Z=Z,  linpar=linpar, Psi=Psi)
}

ffa_gen_L <- function(p, q) {
  matrix(rnorm(p*q), p, q)
}

ffa_gen_Psi <- function(p) {
  runif(p, .5, 2)
}

#' Computes the maximum likelihood estimator for factor analysis
#' @export
ffa <- function(Y, q, maxiter=100, eps=1e-4, savepath=F, verbose=T, iteratively_update_Psi = T, rotate_updates = F){
  stopifnot(is.matrix(Y))

  Miss <- ffa_get_Miss(Y)
  Y <- scale(Y, scale=F)
  p <- ncol(Y)

  beta <- attr(Y, "scaled:center")
  Y_vars <- ffa_comp_Y_vars(Y, Miss)
  A <- diag(1, p, q)
  Psi <- rep(1, p)
  if (savepath) path <- list()

  if (!is.null(Miss)) {
    Y[Miss] <- 0
  }


  for (i in 1:maxiter) {
    Psi.old <- Psi
    A.old <- A
    Zdat <- ffa_est_Z(Y, A, Psi, Miss)
    Z <- Zdat$Z
    A <- ffa_est_A(Y, Z, covZ=Zdat$covZ, covZ.neg = Zdat$covZ.neg, Miss)

    if(iteratively_update_Psi) Psi <- ffa_est_Psi(Y_vars, A, Miss, Z)

    if (!is.null(Miss)) {
      # TODO: this is extremely inefficient (computing the whole matrix just for a few missing value)... do this better
      Y[Miss] <- (Z %*% t(A))[Miss]
    }

    if (savepath) {
      path[[i]] <- list(A = as.vector(A), Psi = Psi)
    }

    if (iteratively_update_Psi) {
      crit <- sum(Psi - Psi.old)^2 / sum(Psi^2 + Psi.old^2)
    } else {
      crit <- eps + 1
    }
    if ((i > 4) && crit < eps) break()
    if (rotate_updates) A <- varimax(A)$loadings
  }
  if(!iteratively_update_Psi) Psi <- ffa_est_Psi(Y_vars, A, Miss, Z)

  if (savepath) {
    path <- sapply(c("A", "Psi"), function(est) {
      do.call(rbind, sapply(seq_along(path), function(i) path[[i]][[est]], simplify=F))
    })
  }
  list(A=A, Psi=Psi, beta=beta, Z=Z, niter=i, path=if(savepath)path)
}

ffa_est_Z <- function(Y, A, Psi, Miss=NULL){
  # For now we estimate missing values by iterative imputation, but this should be better:
  Miss <- NULL
  if (!is.null(Miss)) {
    Y[Miss] <- 0
    # this is the same as ignoring the missing values for each Y[i,]:
    # K <- ffa_comp_K(A, Psi)
    # Z <- t(sapply(1:n, function(i) Y[i,!Miss[i,]] %*% K[!Miss[i,],]))
  }
  K <- ffa_comp_K(A, Psi)
  Z <- Y %*% K

  if (!is.null(Miss)) {
    Z <- (Y / (rowMeans(!Miss))) %*% K
    Z <- scale(Z, scale=F) # TODO: think of this...
  }

  SZZ <- t(Z) %*% Z/nrow(Z)
  CSZZ <- chol(SZZ)

  # compute the target covariance for Z
  covZ.neg <- t(A) %*% K
  covZ <- solve(covZ.neg)
  CcovZ <- chol(covZ)

  # rescale
  # TODO: test: ALL EQUAL, ZT Z / nrow(Z), covZ
  Q <- solve(CSZZ) %*% CcovZ

  Z <- Z %*% Q
  list(Z=Z, covZ = covZ, covZ.neg = covZ.neg)
}

ffa_est_A <- function(Y, Z, covZ=NULL, covZ.neg = NULL, Miss=NULL){
  # For now we estimate missing values by iterative imputation, but this should be better:
  Miss <- NULL

  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(Z)
  if (!is.null(Miss)) {
    Y[Miss] <- 0
  }

  if (is.null(covZ)) {
    covZ <- (t(Z) %*% Z)/nrow(Z)
  }

  if (is.null(covZ.neg)) {
    covZ.neg <- solve(covZ)
  }

  if (!is.null(Miss)) {
    # The following is the same as:
    # A <- t(sapply(1:p, function(j){
    #   t(Y[!Miss[,j],j, drop=F]) %*% Z[!Miss[,j],] %*% (covZ.neg / sum(!Miss[,j]))
    # }))

    if (n < p) {
      A <- ((t(Y)/colSums(!Miss)) %*% Z) %*% (covZ.neg)
    } else {
      A <- (t(Y)/colSums(!Miss)) %*% (Z %*% (covZ.neg))
    }

  } else {
    if (n < p) {
      A <- (t(Y) %*% Z) %*% (covZ.neg / nrow(Y))
    } else {
      A <- t(Y) %*% (Z %*% (covZ.neg / nrow(Y)))
    }
  }
  A
}


ffa_est_Psi <- function(Y_vars, A, Miss, Z, Y) {
  # For now we estimate missing values by iterative imputation, but this should be better:
  Miss <- NULL

  if(is.null(Miss)) {
    Psi <- Y_vars - rowSums(A^2) #*colMeans(!Miss)
  } else {
    Psi <- Y_vars - rowSums(A^2)
  }
  too_small <- Psi < 1e-1
  if(any(too_small)){
    # warning("Psi too small: takes value of 1e-2 instead.")
    Psi[too_small] <- 1e-1
  }
  Psi
}

ffa_get_Miss <- function(Y) {
  isna <- is.na(Y)
  if (any(isna)) {
    if (any(rowSums(isna) == ncol(Y))) stop("Some rows have no non-NA values.")
    if (any(colSums(isna) == nrow(Y))) stop("Some columns have no non-NA values.")
    isna
  } else {
    NULL
  }
}

ffa_comp_Y_vars <- function (Y, Miss) {
  if (!is.null(Miss)){
    colMeans(Y^2, na.rm=T)
  } else {
    colMeans(Y^2)
  }
}

ffa_comp_K <- function(A, Psi) {
  q <- ncol(A)
  Ap <- A/Psi
  ((Ap) %*% solve(t(A) %*% Ap + diag(q)))
}

ffa_comp_crit <- function(a, b){
  sum(abs(a-b))/(sum(abs(a)) + sum(abs(b)))
}

ffa_error <- function(A, target, rotate=F){
  if(rotate) A <- psych::Procrustes(A, target)$loadings
  norm(A-target, "F")/(norm(A, "F") + norm(target, "F"))
}

if(0){
  n <- 1000
  p <- 100
  q <- 10

  dat <- ffa_gen_Y(n, p, q)

  Y <- dat$Y
  Y[runif(prod(dim(Y)))<.5] <- NA
  fit <- ffa(Y, q, savepath = T, verbose=T, maxiter=5, eps=1e-5, iteratively_update_Psi = 2)

  ts.plot(fit$path$A[,1:100])

  plot(dat$A[1:1000], psych::Procrustes(fit$A, dat$A)$loadings[1:1000]); abline(0,1,col=2)

  ffa_ffa_error(dat$A, fit$A, rotate = T)
}
