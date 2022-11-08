require(Matrix)

##' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}


initialize_hessian_AB <- function(fg) {
  fg$parameters$hessian <- lapply(1:fg$dimensions$p, function(na)diag(10, fg$dimensions$q + fg$dimensions$k))
}


update_hessian_AB <- function(hessian_old, hessian_new, weight) {
  warning("DO NOT UPDATE HESSIAN; DO NOT INCLUDE TIME DEPENDENCE THAT WAY")
  lapply(seq_along(hessian_old), function(i){
    weight * hessian_old[[i]] + (1-weight) * hessian_new[[i]]
  })
}


mult_invHessian_dAB <- function(dAB, hessian_AB) {
  prod <- sapply(seq_along(hessian_AB), function(j) {
    as.vector(solve(hessian_AB[[j]], dAB[j,]))
  }, simplify=F)
  do.call(rbind, prod)
}


# Compute the hessian for AB
simulate_hessian_AB <- function(fg){
  Z <- scale(gen_Z(fg$dimensions$n, fg$dimensions$q), scale=F, center=T)
  warning("Z has been rescaled in simulate_hessian_AB")
  linpar <- compute_linpar(Z, fg$parameters$A, fg$X, fg$parameters$B)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, fg$families)
  ZX <- cbind(Z, fg$X)
  sapply(1:fg$dimensions$p, function(j) {
    -(t(ZX) %*% (ZX*(linpar_bprimeprime[,j])))/(fg$parameters$phi[j]*nrow(ZX))
  }, simplify=F)
}
