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
  H_sam   <- compute_psi_AB_hessian(ZX_join(Z, X), phi=parameters$phi, linpar_bprimeprime = linpar_bprimeprime, Miss = Miss)
  # H_sam   <- lapply(H_sam, function(H) diag(diag(H)))



  # AB <- (psi_sim$AB - psi_sam$AB)/nrow(Y)
  # compute the updates
  AB <- compute_hessian_x_psi(psi_sam$AB - psi_sim$AB, H_sam)
    # The above is equivalent to
  AB <- AB_separate(AB, ncol(parameters$Z))

  psi_update <- list(
    A = AB$A,
    # B = parameters_sim$B - parameters_sam$B + AB$B,
    B = AB$B,
    phi = psi_sim$phi - psi_sam$phi,
    # phi = parameters$phi - psi_sam$phi,
    Z = parameters$Z - parameters_sam$Z
  )

  psi_update
}

# returns the (simple) psi functions for A, B, phi
compute_psi_simple <- function(Y, X, Z, parameters, families, Miss, compute_hessian = T) {

  # overhead computations
  if(!is.null(parameters$B)){
    XB <- X %*% t(parameters$B)
  } else {
    XB <- NULL
  }

  # linpar <- compute_linpar(Z, parameters$A, XB=XB)
  # linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  # linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)

  #compute psi_star_AB
  psi_AB  <- compute_psi_simple_AB(Y, Z, X)
  psi_phi <- compute_psi_simple_phi(Y, parameters, families)

  # if(compute_hessian) {
    # psi_AB_hessian <- compute_psi_AB_hessian(ZX, parameters$phi, linpar_bprimeprime, Miss)
  # } else {
    # psi_AB_hessian <- NULL
  # }
  c(psi_AB, list(phi = psi_phi))
}

compute_psi_simple_AB <- function(Y, Z, X) {
  ZX <- ZX_join(Z, X)
   # AB <- log(t(Y+.1)) %*% ZX
  # AB <- t(Y) %*% ZX
  # Y <- scale(log(Y+.1), scale=F)
  AB <- t(Y) %*% ZX
  # Y <- scale(Y, scale=F)
  # A <- t(Y) %*% Z
  # B <- t(matrix(attr(Y, "scaled:center"), ncol=ncol(Y), nrow=nrow(Y), byrow=T)) %*% X # completely inefficient
  # AB <- cbind(A, B)
  # AB <- t(scale(Y, scale=F)) %*% ZX
  # warning("RESCALING HERE LOSES INFORMATION FOR BETA -- THIS IS OBVIOUS FOR GAUSSIAN VARIABLE")

  list(AB=AB)
}

compute_psi_simple_phi <- function(Y, parameters, families) {
  psi_phi <- parameters$phi
  if(length(families$id$gaussian) > 0) {
    # wrong psi_phi[families$id$gaussian] <- colMeans(Y[,families$id$gaussian, drop=F]^2) - rowSums(parameters$A[families$id$gaussian,, drop=F]^2)
    psi_phi[families$id$gaussian] <- colMeans(scale(Y[,families$id$gaussian, drop=F], scale=F)^2)# - rowSums(parameters$A[families$id$gaussian,, drop=F]^2)
    # psi_phi[families$id$gaussian] <- colMeans(Y*linpar) # - colMeans(linpar**2)
  }
  psi_phi/10
}

compute_psi_hessian <- function(Z, linpar_bprimeprime, phi){
  sapply(1:length(phi), function(j) {
    -(t(Z) %*% (Z*(linpar_bprimeprime[,j])))/phi[j]
  }, simplify=F)
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
  psi <- compute_psi_simple(fg$Y, fg$X, fg$parameters$Z, fg$parameters, fg$families, fg$Miss)
  gradient <- compute_gradients_simple(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss)

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
  # Testing behavior of rescale
  devtools::load_all()
  set.seed(1234)
  poisson  <-0
  gaussian <- 10
  binomial <- 10
  q <- 1
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(120303)
  fg <- gen_fastgllvm(nobs=1000, p=p, q=q, family=family, phi=3*(1:p)/p, k=1, intercept=T, miss.prob = 0, scale=1)

  param1 <- fg$parameters
  param1$Z <- with(param1, compute_Z(fg$Y, fg$X, fg$parameters, start=fg$parameters$Z)$Z)

  param2 <- fg$parameters
  param2$B <- param2$B +
  param2$Z <- with(param2, compute_Z(fg$Y, fg$X, fg$parameteres, start=fg$parameters$Z)$Z)

  plot(param1$Z, param2$Z); abline(0,1)

  param1 <- recenter(param1, 1)
  param2 <- recenter(param2, 1)

  plot(param1$Z, param2$Z); abline(0,1)

  param1$B
  param2$B
  all.equal(param1$B, param2$B)



}
