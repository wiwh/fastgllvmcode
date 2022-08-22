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

compute_gradients_simple <- function(Y, X, parameters, families, Miss, ...) {
  parameters_old <- parameters

  if(!is.null(parameters$B) && all(X[,1]==1)) rescale.B=1 else rescale.B=FALSE

  # parameters <- rescale(parameters, rescale.A=T, rescale.B=rescale.B, target.cov = parameters$covZ)

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


  Z0 <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=parameters_old$Z, Miss=Miss)$Zstar
  Z0 <- scale(Z0, scale=F)
  # Z0 <- rescale(parameters)$Z
  # Obtain Zh
  Zh <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar
  Zh <- scale(Zh, scale=F)
  covZ <- cov(Zh)
  # Zh <- rescale(parameters)$Z

  # compute psi on this
  psi_sam <- compute_psi_star_known_Z(Y, X, Z0, parameters, families, Miss, compute_hessian=T)
  psi_sim <- compute_psi_star_known_Z(Y_sim$Y, X, Zh, parameters, families, Miss, compute_hessian=T)


  # Compute both updates
  AB_update_sam <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)
  AB_update_sim <- compute_hessian_x_psi(psi_sim$psi_AB, psi_sim$psi_AB_hessian)
  AB_update <- AB_update_sam - AB_update_sim

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  A_update   <- parameters_old$A - parameters$A + AB_update$A
  B_update   <- parameters_old$B - parameters$B + AB_update$B
  phi_update <- psi_sim$psi_phi  - psi_sam$psi_phi
  Z_update   <- parameters_old$Z - psi_sam$Z
  covZ_update <- parameters$covZ - covZ

  list(A = A_update, B= B_update, phi=phi_update, Z=Z_update, covZ=covZ_update)
}
