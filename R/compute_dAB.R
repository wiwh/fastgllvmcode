delete_initialize_parameters_simple <- function(parameters, dimensions) {
  with(dimensions, {
    if (is.null(parameters$A)) parameters$A <- matrix(0, p, q)
    if (is.null(parameters$B)) parameters$B <- matrix(0, p, k)
    if (is.null(parameters$phi)) parameters$phi <- rep(1, p)
    # if (is.null(Z)) Z <- matrix(0, n, q)
    if (is.null(parameters$covZ)) parameters$covZ <- diag(q)

    # Control parameters
    stopifnot(is.matrix(parameters$A) && dim(parameters$A) == c(p, q))
    stopifnot(is.matrix(parameters$B) && dim(parameters$B) == c(p, k))
    stopifnot(is.vector(parameters$phi) && length(parameters$phi) == p)
    # stopifnot(is.matrix(Z) && dim(Z) == c(n, q))
    stopifnot(is.matrix(parameters$covZ) && dim(parameters$covZ) == c(q,q))
    parameters
  })
}

initialize_gradients_simple <- function(parameters) {
  sapply(parameters, function(par) par*0, simplify=F)
}

compute_dAB_centered <- function(fg, controls, hessian) {
  # Update main values and compute gradient of the sample
  fg$Z <- compute_Z(fg, start=fg$Z, maxit=10)$Z # TODO: maybe reduce to a maxit of 1?
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  fg <- compute_mean(fg)

  # experimental vvvvv
  if (controls$method=="approx") {
    fg_simulated <- fg
    fg_simulated$Z <- scale(gen_Z(fg$dim$n, fg$dimensions$q), scale=F)
    fg_simulated <- compute_mean(fg_simulated)

    dAB <- t(fg_simulated$mean) %*% cbind(fg_simulated$Z, fg_simulated$X) -t(fg$Y) %*% cbind(fg$Z, fg$X)
    dAB <- AB_separate(dAB, fg$dimensions)
    # these are placeholders
    dphi <- rep(0, fg$dimensions$p)
    covZ <- diag(1, fg$dimensions$q)

  # experimental ^^^^^
  # experimental vvvvv
  } else if (controls$method == "rescaled") {
    # Simulate a sample
    fg_simulated <- simulate(fg, return_fastgllvm=T)
    fg_simulated$Z <- compute_Z(fg_simulated, start=fg_simulated$Z, maxit=10)$Z

    # fg_simulated_rescaled <- rescale(fg_simulated, rescale.A = FALSE, rescale.B = FALSE, target.cov = )
    fg_simulated_rescaled <- fg_simulated
    fg_simulated_rescaled <- compute_mean(fg_simulated_rescaled)

    # Compute the centered gradients on both the sample and simulated fastgllvm objects
      fg_rescaled <- rescale(fg, rescale.A= FALSE, rescale.B= TRUE, target.cov = cov(fg_simulated$Z))
    fg_rescaled <- compute_mean(fg_rescaled)

    dAB_sample <- compute_dAB(fg_rescaled, controls$method)
    dAB_simulated <- compute_dAB(fg_simulated_rescaled, controls$method)


    if (is.null(hessian)) {
      dAB <- dAB_simulated - dAB_sample # multiplied by -1 because the hessian is absent and we want to maximize, not minimize
    } else {
      dAB <- mult_invHessian_dAB(dAB_sample - dAB_simulated, hessian)
    }

    if (controls$use_signs) {
      dAB <- sign(dAB)
    }
    dAB <- AB_separate(dAB, fg$dimensions)

    dphi_sample <- compute_phi(fg_rescaled)
    dphi_simulated <- compute_phi(fg_simulated_rescaled)

    dphi <- dphi_simulated - dphi_sample

    covZ <- cov(fg_simulated$Z)

  # experimental ^^^^^
  } else {
    # Simulate a sample
    fg_simulated <- simulate(fg, return_fastgllvm=T)
    fg_simulated$Z <- compute_Z(fg_simulated, start=fg_simulated$Z, maxit=10)$Z
    fg_simulated <- compute_mean(fg_simulated)

    # Compute the centered gradients on both the sample and simulated fastgllvm objects
    dAB_sample <- compute_dAB(fg, controls$method)
    dAB_simulated <- compute_dAB(fg_simulated, controls$method)


    if (is.null(hessian)) {
      dAB <- dAB_simulated - dAB_sample # multiplied by -1 because the hessian is absent and we want to maximize, not minimize
    } else {
      dAB <- mult_invHessian_dAB(dAB_sample - dAB_simulated, hessian)
    }

    if (controls$use_signs) {
      dAB <- sign(dAB)
    }
    dAB <- AB_separate(dAB, fg$dimensions)

    dphi_sample <- compute_phi(fg)
    dphi_simulated <- compute_phi(fg_simulated)

    dphi <- dphi_simulated - dphi_sample

    covZ <- cov(fg_simulated$Z)
  }


  # dcovZ <- cov(fg$Z) - cov(fg_simulated$Z)


  # TODO: compare to the old method below, they need to get the same thing!!
  list(dA = dAB$A, dB=dAB$B, fg=fg, covZ=covZ, dphi=dphi)#, dcovZ=dcovZ)
}


compute_dAB <- function(fg, method) {
  with(fg, {
    if (dimensions$k == 0) {
      ZX <- Z
    } else {
      ZX <- cbind(Z, X)
    }

    if ( method == "simple" ) {
      dAB <- t(Y) %*% (ZX/dimensions$n)
    } else if (method == "full") {
      dAB <- t(Y - mean) %*% (ZX/dimensions$n)
    } else if (method == "approx") {
      dAB <- t(mean) %*% ZX /dimensions$n

    # TODO: remove the method
    } else if (method == "rescaled") {
      dAB <- t(Y - mean) %*% (ZX/dimensions$n)
    } else {
      stop ("Unkown method `method`")
    }
    dAB
  })
}

compute_gradients_simple_old <- function(Y, X, parameters, families,  ...) {
  parameters_sim <- parameters_sam <- parameters

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

  # Compute psi for sam
  parameters_sam$Z <- scale(compute_Z(Y, X, parameters_sam, families, start=parameters_sam$Z)$Z, scale=F, center=T) # DO NOT RESCALE SCALE, CENTER IS OK
  # parameters_sam$Z <- rescale(parameters_sam, T, 1)$Z
  psi_sam <- compute_psi_simple(Y, X, parameters_sam$Z, parameters_sam, families)

  # Compute psi for sim
  parameters_sim$Z <- scale(compute_Z(Y, X, parameters_sim, families, start=Y_sim$Z)$Z, scale=F, center=T) # DO NOT RESCALE SCALE, CENTER IS OK
  # parameters_sim$Z <- rescale(parameters_sim, T, 1)$Z
  psi_sim <- compute_psi_simple(Y_sim$Y, X, parameters_sim$Z, parameters_sim, families)

  # compute independenz Z for hessian:
  Z <- scale(gen_Z(nrow(Y), ncol(parameters$A)), scale=F) # this prevents bias for computing the hessian: it is independent of the rest.
  linpar <- compute_linpar(Z, parameters$A, X, parameters$B)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)
  H_sam   <- compute_psi_AB_hessian(ZX_join(Z, X), phi=parameters$phi, linpar_bprimeprime = linpar_bprimeprime)
  # H_sam   <- lapply(H_sam, function(H) diag(diag(H)))



  # AB <- (psi_sim$AB - psi_sam$AB)/nrow(Y)
  # compute the updates
  AB <- compute_hessian_x_psi(psi_sam$AB - psi_sim$AB, H_sam)
    # The above is equivalent to
  AB <- AB_separate(AB, ncol(Z))

  psi_update <- list(
    A = AB$A,
    # B = parameters_sim$B - parameters_sam$B + AB$B,
    B = AB$B,
    phi = psi_sim$phi - psi_sam$phi,
    # phi = parameters$phi - psi_sam$phi,
    Z = Z - parameters_sam$Z
  )

  psi_update
}

compute_psi_simple_phi <- function(Y, parameters, families) {
  psi_phi <- parameters$phi * 0
  if(length(families$id$gaussian) > 0) {
    # wrong psi_phi[families$id$gaussian] <- colMeans(Y[,families$id$gaussian, drop=F]^2) - rowSums(parameters$A[families$id$gaussian,, drop=F]^2)
    psi_phi[families$id$gaussian] <- colMeans(scale(Y[,families$id$gaussian, drop=F], scale=F)^2)# - rowSums(parameters$A[families$id$gaussian,, drop=F]^2)
    # psi_phi[families$id$gaussian] <- colMeans(Y*linpar) # - colMeans(linpar**2)
  }
  psi_phi/10
}


if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 0
  gaussian <- 10
  binomial <- 0
  q <- 2
  k <- 1
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(1030)
  fg <- gen_fastgllvm(nobs=100, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, intercept=T, miss.prob = 0, scale=1)
  psi <- compute_psi_simple(fg$Y, fg$X, fg$Z, fg$parameters, fg$families)
  gradient <- compute_gradients_simple(fg$Y, fg$X, fg$parameters, fg$families)

  # this must have expectaiton 0 under the true model. check this
  sim <- sapply(1:1000, function(i){
    set.seed(i)
    if(i%%10==0 || i==1) cat("\n", i)
    fg <- gen_fastgllvm(nobs=100, p=p, q=q, A= fg$parameters$A, B = fg$parameters$B, phi=fg$parameters$phi, family=family, k=1, intercept=T, miss.prob = 0, scale=1)
    psi <- compute_gradients_simple(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss)
    psi
  }, simplify=F)

  sim <- sapply(names(sim[[1]]),
    function(parname) do.call(rbind, sapply(sim, function(simi) as.vector(simi[[parname]]),
    simplify=F)), simplify=F)

  image(cov(cbind(sim$A, sim$B))) # super big

  boxplot(sim$A, outline=F); abline(h=0, col=2)
  boxplot(sim$B,outline=F); abline(h=0, col=2)
  points(colMeans(sim$B), col=2)
  boxplot(sim$phi);abline(h=0, col=2)
}

if(0) {
  # Testing behavior of the different scores
  devtools::load_all()
  set.seed(1234)
  poisson  <- 100
  gaussian <- 0
  binomial <- 0
  q <- 1
  k <- 1
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(12395)
  fg <- gen_fastgllvm(nobs=100, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, intercept=T, miss.prob = 0, scale=1)

  fg$hessian <- simulate_hessian_AB(fg)
  for(i in 1:100){
    fg$hessian <- update_hessian_AB(fg$hessian, simulate_hessian_AB(fg), weight=.95)
  }
  # fg$hessian <- lapply(hessian, function(na) diag(1, ncol(na)))

  method <- "full"
  use_signs <- T
  hessian <- NULL

  set.seed(123)
  sims <- t(sapply(1:1000, function(na) {
    fg <- simulate(fg, return_fastgllvm = TRUE)
    controls <- list(method=method, use_signs=use_signs)
    as.vector(unlist(compute_dAB_centered(fg, controls=controls, hessian=hessian)[c("dA", "dB")]))
  }))
  {
  par(mfrow=c(2,1))
  boxplot(sims, outline=T, main=sqrt(mean(colMeans(sims)^2)))
  points(colMeans(sims), col=2)
  abline(h=0,col=2, lty=2)
  boxplot(sims, outline=F)
  points(colMeans(sims), col=2)
  abline(h=0,col=2, lty=2)
  par(mfrow=c(1,1))
  }
}
