initialize_parameters_simple <- function(parameters, dimensions) {
  with(dimensions, {
    if (is.null(parameters$A)) parameters$A <- matrix(0, p, q)
    if (is.null(parameters$B)) parameters$B <- matrix(0, p, k)
    if (is.null(parameters$phi)) parameters$phi <- rep(1, p)
    if (is.null(parameters$Z)) parameters$Z <- matrix(0, n, q)
    if (is.null(parameters$covZ)) parameters$covZ <- diag(q)

    # Control parameters
    stopifnot(is.matrix(parameters$A) && dim(parameters$A) == c(p, q))
    stopifnot(is.matrix(parameters$B) && dim(parameters$B) == c(p, k))
    stopifnot(is.vector(parameters$phi) && length(parameters$phi) == p)
    stopifnot(is.matrix(parameters$Z) && dim(parameters$Z) == c(n, q))
    stopifnot(is.matrix(parameters$covZ) && dim(parameters$covZ) == c(q,q))
    parameters
  })
}

initialize_gradients_simple <- function(parameters) {
  sapply(parameters, function(par) par*0, simplify=F)
}

compute_gradients_simple <- function(Y, X, parameters, families, Miss, debiase) {
  A_old <- parameters$A
  Z_old <- parameters$Z
  # begin by rescaling
  resc <- rescale(parameters$Z, parameters$A, target.cov=parameters$covZ)
  parameters$A <- resc$A
  parameters$Z <- resc$Z

  # update_A <- A_old - parameters$A
  # Compute zhat on sample
  Z0 <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=parameters$Z, Miss=Miss)$Zstar
  # rescale
  # A_old <- parameters$A
  # TODO: optimize the rescaling for B too...
  if(!is.null(parameters$B) && all(X[,1]==1)) {
    Z0 <- scale(Z0, scale=F)
    parameters$B[,1] <- parameters$B[,1]  - as.vector(parameters$A %*% attr(Z0, "scaled:center"))
  }
  resc <- rescale(Z0, parameters$A, target.cov=parameters$covZ)
  parameters$A <- resc$A
  Z0 <- resc$Z

  # Generate sim
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y)
  )

  # Obtain Zh
  Zh <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar
  # update covZ
  # covZ <- .9*parameters$covZ + .1 * cov(Zh)
  covZ <- cov(Zh)
  covZ_update <- parameters$covZ - covZ # this will be substracted from parameters$covZ


  # rescale Zh
  resc <- rescale(Zh, target.cov = parameters$covZ)
  Zh <- resc$Z

  # compute psi on this
  psi_sam <- compute_psi_star_known_Z(Y, X, Z0, parameters, families, Miss, compute_hessian=T)
  psi_sim <- compute_psi_star_known_Z(Y_sim$Y, X, Zh, parameters, families, Miss, compute_hessian=T)


  # update
  if(debiase) {
    # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)
    AB_update_sam <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)
    AB_update_sim <- compute_hessian_x_psi(psi_sim$psi_AB, psi_sim$psi_AB_hessian)
    AB_update <- AB_update_sam - AB_update_sim
  } else {
    AB_update <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)
  }
  AB_update <- AB_separate(AB_update, ncol(parameters$A))
  phi_update <- psi_sim$psi_phi - psi_sam$psi_phi

  Z_update <- Z_old - psi_sam$Z # TODO: or take from the parameter update

  list(A = AB_update$A + A_old - parameters$A, B= AB_update$B, phi=phi_update, Z=Z_update, covZ=covZ_update)
}
