# this is to compute dz / dlambda for inference purposes.

compute_psiz_da <- function(A, phi, linpar, families, Miss) {
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar, families)

  p <- nrow(A)

  if(!is.null(Miss)) {
    # TODO: do this in parent!
    linpar_bprimeprime[Miss] <- 0
  }

  lapply(1:ncol(linpar_bprimeprime), function(j) {

  })

  # TODO: this should be done with tensor products
  if(p > 1) {
    lapply(1:ncol(linpar_bprimeprime), function(j) {
      - (t(A) %*% (A * linpar_bprimeprime[j,] / phi) + diag(q))
    })
  } else {
    - (t(t(linpar_bprimeprime)*(as.vector(A)/phi)) %*% A + 1)  # TODO: explain this vector trick
  }
}


if(0) {
  devtools::load_all()
}
