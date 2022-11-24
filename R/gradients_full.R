initialize_parameters_full <- function(parameters, dimensions) {
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

initialize_gradients_full <- function(parameters) {
  sapply(parameters, function(par) par*0, simplify=F)
}

compute_gradients_full <- function(Y, X, parameters, families, Miss, debiase) {
  A_old <- parameters$A
  Z_old <- parameters$Z
  # begin by rescaling
  resc <- rescale.old(parameters$Z, parameters$A, target.cov=parameters$covZ)
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
  resc <- rescale.old(Z0, parameters$A, target.cov=parameters$covZ)
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
  # TODO: generate without miss ! or with miss?
  # Obtain Zh
  Zh <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar
  # update covZ
  # covZ <- .9*parameters$covZ + .1 * cov(Zh)
  covZ <- cov(Zh)
  covZ_update <- parameters$covZ - covZ # this will be substracted from parameters$covZ


  # rescale.old Zh
  resc <- rescale.old(Zh, target.cov = parameters$covZ)
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
  phi_update <- (psi_sim$psi_phi - psi_sam$psi_phi)

  Z_update <- Z_old - psi_sam$Z # TODO: or take from the parameter update

  list(A = AB_update$A + A_old - parameters$A, B= AB_update$B, phi=phi_update, Z=Z_update, covZ=covZ_update)
}

compute_gradients_full <- function(Y, X, parameters, families, Miss, debiase) {
  parameters_sam <- parameters_sim <- parameters
  parameters_sam$Z <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=parameters$Z, Miss=Miss)$Zstar

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
  parameters_sim$Z <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar

  # compute psi on this
  psi_sam <- compute_psi_star_known_Z(Y, X, parameters_sam$Z, parameters, families, Miss, compute_hessian=T)
  psi_sim <- compute_psi_star_known_Z(Y_sim$Y, X, parameters_sim$Z, parameters, families, Miss, compute_hessian=T)


  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)
  AB_update_sam <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)
  AB_update_sim <- compute_hessian_x_psi(psi_sim$psi_AB, psi_sim$psi_AB_hessian)
  AB_update <- AB_update_sam - AB_update_sim

  AB_update <- AB_separate(AB_update, ncol(parameters$A))
  phi_update <- (psi_sim$psi_phi - psi_sam$psi_phi)

  Z_update <- parameters$Z - parameters_sam$Z


  list(A = AB_update$A, B= AB_update$B, phi=phi_update, Z=Z_update)
}




#  used only for testing, AB is never computed at each iteration...
compute_AB_update <- function (Y, Z, X, B, A, phi, families, Miss=NULL) {
  linpar <- compute_linpar(Z, A, X=X, B=B)
  linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families) # should not appear here

  ZX <- ZX_join(Z, X)
  psi_AB <- compute_psi_AB(Y, ZX, phi, linpar_bprime, Miss=Miss)
  psi_AB_hess <- compute_psi_AB_hessian(ZX, phi = phi, linpar_bprimeprime = linpar_bprimeprime, Miss=Miss)

  AB <- cbind(A,B)
  AB <- AB - compute_hessian_x_psi(psi_AB, psi_AB_hess)

  AB
}

#This returns a list of high-level updates for the parameter
compute_psi_star_known_Z <- function(Y, X, Z, parameters, families, Miss, compute_hessian, Z.maxit=100) {
  # overhead computations
  if(!is.null(parameters$B)){
    XB <- X %*% t(parameters$B)
  } else {
    XB <- NULL
  }

  linpar <- compute_linpar(Z, parameters$A, XB=XB)
  linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)

  ZX <- ZX_join(Z, X)
  #compute psi_star_AB
  psi_AB <- compute_psi_AB(Y, ZX, parameters$phi, linpar_bprime, Miss=Miss)

  psi_phi <- compute_psi_phi(Y, parameters$phi, linpar_bprime, families, Miss)

  if(compute_hessian) {
    psi_AB_hessian <- compute_psi_AB_hessian(ZX, parameters$phi, linpar_bprimeprime, Miss)
  } else {
    psi_AB_hessian <- NULL
  }


  list(psi_AB = psi_AB, psi_phi = psi_phi, psi_AB_hessian = psi_AB_hessian, linpar=linpar, Z=Z)
}


compute_hessian_x_psi <- function(AB_psi, AB_hessian) {
  prod <- sapply(seq_along(AB_hessian), function(j) {
    as.vector(solve(AB_hessian[[j]], AB_psi[j,]))
  }, simplify=F)
  do.call(rbind, prod)
}

compute_psi_AB <- function(Y, ZX, phi, linpar_bprime, Miss){
  if (!is.null(Miss)) {
    Y[Miss] <- 0 # this is correct... check the NA rmd file for examples
    linpar_bprime[Miss] <- 0
  }
  (t(Y - linpar_bprime)/phi) %*% ZX
}


compute_psi_AB_hessian <- function (ZX, phi, linpar_bprimeprime, Miss) {
  if (!is.null(Miss)) {
    linpar_bprimeprime[Miss] <- 0 # do at the parent lvl... this trick is OK, check NA.RMd
  }
  # compute a list of all hessians
  # TODO if Z and X are standardized, this matrix is diagonal.
  sapply(1:length(phi), function(j) {
    -(t(ZX) %*% (ZX*(linpar_bprimeprime[,j])))/phi[j]
  }, simplify=F)
}




# This is only for testing purposes, supposed to be equal to the above
compute_psi_AB_test <- function(Y, ZX, phi, linpar_bprime, Miss){
  if(!is.null(Miss)) {
    Y[miss] <- 0
  }
  diag(1/phi) %*% t(Y - linpar_bprime) %*% ZX
}

compute_psi_phi <- function (Y, phi, linpar_bprime, families, Miss) {
  psi_phi <- rep(0, length(phi))
  # update gaussian
  if(length(families$id$gaussian > 0)) {
    id <- families$id$gaussian
    if(is.null(Miss)) {
      # psi_phi[id] <- (colMeans((Y[,id] - linpar_bprime[,id])**2) - phi[id]) / (2 * phi[id]**2) # this is theoretically correct, but badly behaved
      psi_phi[id] <- (colMeans(scale((Y[,id] - linpar_bprime[,id]), scale=F)**2) - phi[id]) / 5  # this is rescaled appropriately..
      # psi_phi[id] <- t(Y[,id]) %*% linpar_bprime[,id]/nrow(Y[,id])
    } else {
      Y[,id][Miss[,id]] <- 0
      # psi_phi[id] <- (colSums((Y[,id] - linpar_bprime[,id])**2) - phi[id]) / (2 * phi[id]**2) # this is theoretically correct, but badly behaved
      psi_phi[id] <- (colSums(scale((Y[,id] - linpar_bprime[,id]), scale=F)**2) - phi[id]) / 5 # this is rescaled appropriately
      # rescale
      psi_phi[id] <- psi_phi[id]/(colSums(!Miss[,id]))
    }
  }
  psi_phi
}


compute_A_glm <- function(Y, Z, X, families, maxit=100) {
  ZX <- ZX_join(Z, X)
  fit <- sapply(1:ncol(Y), function(j) {
    glm(Y[,j] ~0 + ZX, family = families$vec[j], control=list(maxit=maxit))$coef
  }, simplify=F)
  do.call(rbind, fit)
}

compute_error <- function(A1, A2, rotate=T) {
  if(any(dim(A1)!=dim(A2))) stop("Dimensions unequal.")
  if(rotate) A1 <- psych::Procrustes(A1, A2)$loadings
  norm(A1 - A2, type="F")/(norm(A1, type="F") + norm(A2, type="F"))
}

# Update experiments
# ---------

if(0) {
  devtools::load_all()
  set.seed(121234)
  poisson  <- 0
  gaussian <- 10
  binomial <- 0
  q <- 2
  p <- poisson + gaussian + binomial
  fg <- gen_fastgllvm(nobs=1000, p=p, q=q, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)),k=1, intercept=T, miss.prob = 0, scale=1)

  #TODO: problem with missing value when q=1
  values <- list(parameters = compute_parameters_initial_values(fg, rescale=F))
  values <- c(values, list(gradients = initialize_gradients(values$parameters)))

  compute_gradients_function <- compute_gradients

  params_hist <- list()
  moving_average <- values$parameters
  for(i in 1:50){
    moving_average_old <- moving_average
    if(i==10) print("Changing to debiase now!")
    if(i < 10) {
      values <- update_parameters(fg$Y, fg$X, values$parameters, values$gradients, fg$families, fg$Miss, alpha=1, beta=0, debiase=T, compute_gradients = compute_gradients_function)
    } else {
      values <- update_parameters(fg$Y, fg$X, values$parameters, values$gradients, fg$families, fg$Miss, alpha=2/(2+i**.8), beta=0.2, debiase=T,  compute_gradients = compute_gradients_function)
    }

    for(k in seq_along(parameters)){
      moving_average[[k]] <- moving_average[[k]] * .9 + values$parameters[[k]] * .1
    }
    err <- compute_error(moving_average$A, moving_average_old$A, rotate=F)
    print(err)
    if(i > 20 & err<1e-3) break()
    # print(values$parameters$covZ)
    if(i%%5 == 1 ) {
      Sys.sleep(1)
      par(mfrow=c(2,1))
      plot(fg$parameters$A, psych::Procrustes(values$parameters$A, fg$parameters$A)$loadings, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=1)
      # plot(fg$Z, values$parameters$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
      plot(fg$parameters$B, values$parameters$B, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2)
      par(mfrow=c(1,1))
    }
    params_hist <- c(params_hist, list(as.vector(values$parameters$A)))
  }
  params_hist <- do.call(rbind, params_hist)
  ts.plot(params_hist)
  points(rep(nrow(params_hist), ncol(params_hist)), as.vector(moving_average$A),col=2)
  #
  # library(gllvm)
  # fit.gllvm <- gllvm(y=fg$Y, x=fg$X, formula=~0+x, num.lv=1, family="binomial", sd.errors=F)
  #
  library(gmf)
  # fit.gmf <- gmf(fg$Y, X=fg$X, family=binomial(), p = q)
  fit.gmf <- gmf(fg$Y, X=fg$X, family=poisson(), p = q)

  library(mirtjml)

  fit.mirtjml <- mirtjml_expr(fg$Y, K= q, tol = 1e-1)
  compute_error(moving_average$A, fg$parameters$A, rotate=T)
  compute_error(fit.gmf$v, fg$parameters$A, rotate=T)
  compute_error(fit.mirtjml$A_hat, fg$parameters$A, rotate=T)

  plot(fg$parameters$A, psych::Procrustes(moving_average$A, fg$parameters$A)$loadings, ylim=c(-5,5), xlim=c(-3,3))
  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=2); abline(0,1,col=1)
  points(fg$parameters$A, psych::Procrustes(fit.mirtjml$A_hat, fg$parameters$A)$loadings, col=3); abline(0,1,col=1)

  plot(fg$parameters$B, moving_average$B)
  points(fg$parameters$B, t(fit.gmf$beta), col=2)
  points(fg$parameters$B, fit.mirtjml$d_hat, col=3)

  compute_error(fg$parameters$B, moving_average$B, rotate=F)
  compute_error(fg$parameters$B, fit.mirtjml$d_hat, rotate=F)
}


# TESTS
# --------
if(0) {

  devtools::load_all()
  set.seed(1423)
  poisson  <- 10
  gaussian <- 10
  binomial <- 10
  p <- poisson + gaussian + binomial
  fg <- gen_fastgllvm(nobs=100, p=p, q=1, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)), k=0, intercept=F, miss.prob = 0)
  zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, Miss=Miss))
  plot(fg$Z, zhat$Zstar)

  # test that the gradient has expectation 0
  sim <- sapply(1:100, function(i){
    fg <- gen_fastgllvm(nobs=100, p=p, q=1, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)), k=0, intercept=F, miss.prob = 0)
    compute_gradients_full(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss, debiase=T)
  })

  sim_A <- do.call(rbind, sapply(sim, function(simi)as.vector(simi$A)))

  AB1 <- with(fg, compute_A_glm(Y, Z, X, families, maxit=1))
  AB2 <- with(fg, compute_A_glm(Y, Z, X, families, maxit=2))

  # Check that starting at AB1 and doing 1 iteration yields AB2
  A <- AB_separate(AB1, fg$dimensions$q, fg$dimensions$k)$A
  B <- AB_separate(AB1, fg$dimensions$q, fg$dimensions$k)$B
  AB <- with(fg, compute_AB_update(Y = Y, Z = Z, X = X, B = B, A=A, phi = parameters$phi, families = families), Miss=Miss)

  # TEST: MUST BE TRUE
  all.equal(AB, AB2)
  # END OF TEST

  # TEST
  linpar_bprime <- with(fg, compute_linpar_bprime(linpar$linpar, families))

  psi_AB1 <- with(fg, compute_psi_AB(Y, cbind(Z,X), parameters$phi, linpar_bprime, Miss=Miss))
  psi_AB2 <- with(fg, compute_psi_AB_test(Y, cbind(Z,X), parameters$phi, linpar_bprime, Miss=Miss))

  all.equal(psi_AB1, psi_AB2)
  # TEST DONE
  # Testing updating parameters... ouchi

}
