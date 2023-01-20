gradient_full <- function (fg) {
  # browser()
  # Compute Z
  fg$Z <- compute_Z(fg, start=fg$Z, maxit=10)$Z # TODO: maybe reduce to a maxit of 1?
  # cat(as.vector(cov(fg$Z)))
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  # Update fg$linpar and fg$mean
  fg <- compute_mean(fg, linpar=NULL, return_object=T)

  dAB <- compute_dAB(fg, method="full")
  dphi <- compute_phi(fg)
  dcovZ <- cov(fg$Z)


  list(AB=dAB, phi=dphi, covZ = dcovZ)
}

gradient_simple <- function(fg) {
  # Update main values and compute gradient of the sample
  fg$Z <- compute_Z(fg, start=fg$Z, maxit=10)$Z # TODO: maybe reduce to a maxit of 1?
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  fg <- compute_mean(fg, return_object=T) # TODO: check if we need this here
  # fg <- rescale(fg, rescale.A=T, rescale.B=T, target.cov=fg$parameters$covZ)

  dAB <- compute_dAB(fg, method="simple")
  dphi <- compute_phi(fg)
  dcovZ <- cov(fg$Z)

  list(AB=dAB, phi=dphi, covZ = dcovZ)
}

