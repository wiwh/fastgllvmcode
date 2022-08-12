# gen_data

ffa_gen_Z <- function(n, q) {
  matrix(rnorm(n*q), n, q)
}

gen_Y <- function(n=NULL, p=p, q=q, A=NULL, Psi=NULL, Z=NULL) {
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
  rep(1, p)
}

#' Computes the maximum likelihood estimator for factor analysis
#' @export
ffa <- function(Y, q, maxiter=100, eps=1e-4, savepath=F, verbose=F){
  stopifnot(is.matrix(Y))

  Y <- scale(Y, scale=F)
  p <- ncol(Y)

  beta <- attr(Y, "scaled:center")
  Y_vars <- colMeans(Y^2)
  A <- diag(1, p, q)
  Psi <- rep(1, p)
  if (savepath) path <- list()

  for (i in 1:maxiter) {
    Psi.old <- Psi
    Zdat <- ffa_est_Z(Y, A, Psi)
    Z <- Zdat$Z
    A <- ffa_est_A(Y, Z, covZ=Zdat$covZ, covZ.neg = Zdat$covZ.neg)

    # in one go
    # browser()
    # K <- ffa_comp_K(A, Psi)
    # KSK <- t(K) %*% (A %*% t(A) + diag(Psi)) %*% K
    # KYK <- t(K) %*% ((t(Y) %*% Y)/nrow(Y)) %*% K
    #
    # A2 <- (t(Y)%*% Y/nrow(Y)) %*% K %*% solve(chol(KYK)) %*% chol(KSK))

    # TEST THAT: all.equal(A1, A2)
    Psi <- ffa_est_Psi(Y_vars, A)
    crit <- ffa_comp_crit(Psi.old, Psi)
    if(verbose)cat("\nIteration: ", i, " - crit: ", crit, ".")

    if (savepath) {
      path[[i]] <- list(A = as.vector(A), Psi = Psi, crit=crit)
    }
    if (crit < eps) break()
  }

  if (savepath) {
    path <- sapply(c("A", "Psi", "crit"), function(est) {
      do.call(rbind, sapply(seq_along(path), function(i) path[[i]][[est]], simplify=F))
    })
  }
  list(A=A, Psi=Psi, beta=beta, Z=Z, niter=i, path=if(savepath)path)
}

ffa_est_Z <- function(Y, A, Psi){
  K <- ffa_comp_K(A, Psi)
  Z <- Y %*% K

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

ffa_est_A <- function(Y, Z, covZ=NULL, covZ.neg = NULL){
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(Z)

  if (is.null(covZ)) {
    covZ <- t(Z) %*% Z
  }

  if (is.null(covZ.neg)) {
    covZ.neg <- solve(covZ)
  }

  if (n < p) {
    A <- (t(Y) %*% Z) %*% (covZ.neg / nrow(Y))
  } else {
    A <- t(Y) %*% (Z %*% (covZ.neg / nrow(Y)))
  }
  A
}

est_A_miss <- function() {

}

est_Z_miss <- function() {

}

ffa_est_Psi <- function(Y_vars, A) {
  Psi <- Y_vars - rowSums(A^2)
  too_small <- Psi < 1e-2
  if(any(too_small)){
    warning("Psi too small: takes value of 1e-2 instead.")
    Psi[too_small] <- 1e-2
  }
  Psi
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
  n <- 100
  p <- 1000
  q <- 10

  dat <- gen_Y(n, p, q)

  fit <- ffa(dat$Y, q, savepath = T, verbose=T, maxiter=100, eps=1e-5)

  ts.plot(fit$path$A[,1:100])
  ts.plot(fit$path$crit)
  ts.plot(fit2$path$A[,1:100])
  ts.plot(fit$path$Psi[,1:100])

  plot(dat$A[1:1000], psych::Procrustes(fit$A, dat$A)$loadings[1:1000]); abline(0,1,col=2)

  ffa_error(dat$A, fit$A, rotate = T)
}
