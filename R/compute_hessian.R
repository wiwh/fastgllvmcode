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
