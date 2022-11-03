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

compute_dAB_centered <- function(fg, method) {
  fg_simulated <- fg_old <- fg

  # Update main values and compute gradient of the sample
  fg$Z <- compute_Z(fg, maxit=10)$Z # TODO: maybe reduce to 1?
  fg$linpar <- with(fg, compute_linpar(Z, parameters$A , X , parameters$B))$linpar
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  fg$mean <- compute_linpar_bprime(fg$linpar, fg$families)

  dAB_sample <- compute_dAB(fg, method)

  # Compute gradient on a simulated sample
  fg_simulated <- simulate(fg, return_fastgllvm=T)
  if (method == "full") fg_simulated$mean <- compute_linpar_bprime(fg_simulated$linpar, fg_simulated$families)
  fg_simulated$Z <- compute_Z(fg_simulated, maxit=10)$Z
  dAB_simulated <- compute_dAB(fg_simulated, method)

  dAB <- mult_invHessian_dAB(dAB_sample - dAB_simulated, fg$hessian)

  # TODO: compare to the old one below, they need to get the same thing!!

  list(dAB = dAB, linpar=fg$linpar, mean=fg$mean)
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
  # TODO: i was here, continue 03.11.2022
  # Testing behavior of rescale
  devtools::load_all()
  set.seed(1234)
  poisson  <- 10
  gaussian <- 10
  binomial <- 10
  q <- 1
  k <- 0
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(120303)
  fg <- gen_fastgllvm(nobs=64, p=p, q=q, family=family, phi=3*(1:p)/p, k=1, intercept=T, miss.prob = 0, scale=1)

  fg$hessian <- compute_hessian_AB(fg)
  for(i in 1:100){
    fg$hessian <- update_hessian_AB(fg$hessian, compute_hessian_AB(fg), weight=.95)
  }
  hessian <- fg$hessian

  set.seed(123)
  sims <- t(sapply(1:1000, function(na) {
    fg <- simulate(fg, return_fastgllvm = TRUE)
    # fg$hessian <- lapply(hessian, function(na) diag(1, ncol(na)))
    fg$hessian <- hessian
    as.vector(compute_dAB_centered(fg, method="full")$dAB)
  }))
  boxplot(sims)
  hist(apply(sims,2,median))
  hist(colMeans(sims))
  apply(sims,2,mean)

  param1 <- fg$parameters
  compute_Z(fg)

  compute_gradients_simple(fg, Z_maxit=1)


  plot(param1$Z, param2$Z); abline(0,1)

  param1 <- recenter(param1, 1)
  param2 <- recenter(param2, 1)

  plot(param1$Z, param2$Z); abline(0,1)

  param1$B
  param2$B
  all.equal(param1$B, param2$B)



}
