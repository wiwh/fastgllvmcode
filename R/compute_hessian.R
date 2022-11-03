# Compute the hessian for AB
compute_hessian_AB <- function(X, dimensions, parameters, families){
  Z <- gen_Z(dimensions$n, dimensions$q)
  linpar <- compute_linpar(Z, parameters$A, X, parameters$B)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)
  ZX <- cbind(Z, X)
  sapply(1:length(phi), function(j) {
    -(t(ZX) %*% (ZX*(linpar_bprimeprime[,j])))/phi[j]
  }, simplify=F)
}
