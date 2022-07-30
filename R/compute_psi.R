#' #' Computes the psi function
#' #'
#' #' @param Y: n*p matrix of responses
#' #' @param Z: n*q matrix of latent variables
#' #' @param X: n*k design matrix
#' #' @param A: p*q loadings matrix
#' #' @param B: p*k coefficients matrix
#' #' @param family: c("gaussian", "negbin", "poisson", "binomial"), or a named
#' #' list of indices. See documentation.
#' #' @param family: known parameters
#' compute_psi <- function (Y, Z, X, A, B, phi, families) {
#' }

# compute_pi <- function(Y, Z, X, A, B, phi, families, maxit=100) {
#   # TODO: this can easily be palalellized
#   sol <- t(sapply(1:ncol(Y), function(j){
#     compute_glm(j, Y, Z, X, A, B, phi, families, maxit=maxit)
#   }))
#   A <- sol[, startsWith(colnames(sol), "A"), drop=F]
#   B <- sol[, startsWith(colnames(sol), "B"), drop=F]
#   phi <- sol[, startsWith(colnames(sol), "phi")]
#
#   list(A=A, B=B, phi=phi)
# }

# compute_glm <- function(j, Y, Z, X, A, B, phi, families, maxit=100) {
#   if(is.null(A) | is.null(B) | is.null(phi)) {
#     fit <- glm(Y[,j]~0+X+Z, family=families$vec[j], maxit=maxit)
#   } else {
#     fit <- glm(Y[,j]~0+X+Z, start=c(B[j,], A[j,]), family=families$vec[j], maxit=maxit)
#   }
#   B_j <- fit$coef[1:ncol(X)]
#   A_j <- fit$coef[(ncol(X)+1):length(fit$coef)]
#   phi_j <- 1
#   if (family=="gaussian") {
#     phi <- var(fit$residuals)
#   }
#   c(A_j=A_j, B_j=B_j, phi_j=phi_j)
# }



# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute psi sample
  # parameters$A[,] <- matrix(rnorm(prod(dim(parameters$A))), nrow(parameters$A), ncol(parameters$A))/5
  # browser()
  psi_sam <- compute_psi_star(Y, X, Z, parameters, families, Miss, compute_hessian=T)
  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )

  # set.seed(304434)
  # sim <- sapply(1:1000, function(aa){
  #   cat("\n", aa)
  #   Y_sim <-  generate_y(
  #       linpar = NULL,
  #       phi = parameters$phi,
  #       families = families,
  #       A = parameters$A,
  #       B = parameters$B,
  #       X = X,
  #       Z = NULL,
  #       nobs = nrow(Y),
  #       Miss = Miss
  #     )
  #   as.vector(compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=F)$psi_AB)
  # }, simplify=F)
  #
  # sim <- do.call(rbind, sim)
  #
  # plot(as.vector(fg$parameters$A), colMeans(sim[,1:1000]))
  # plot(as.vector(fg$parameters$B), colMeans(sim[,1001:2000]))

  # Update AB: each his own hessian
  # psi_sim <- compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=T)
  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian) -
  #              compute_hessian_x_psi(psi_sim$psi_AB, psi_sim$psi_AB_hessian)
  # AB_update <- AB_separate(AB_update, ncol(parameters$A))

  # Update AB: only hessian of Y_0
  psi_sim <- compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=T)
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}


# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters_rescale <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute psi sample
  # parameters$A[,] <- matrix(rnorm(prod(dim(parameters$A))), nrow(parameters$A), ncol(parameters$A))/5
  # browser()
  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )

  # TODO: rescale before or after computing the psi_sim??
  psi_sim <- compute_psi_star_rescaled(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=F, rescale.target=parameters$covZ)
  #update covZ
  psi_sam <- compute_psi_star_rescaled(Y, X, Z, parameters, families, Miss, compute_hessian=T, rescale.target=parameters$covZ)

  parameters$covZ <- 0.9*parameters$covZ + .1 * cov(psi_sim$Z)*(nrow(psi_sim$Z)/(nrow(psi_sim$Z) - 1))

  # Update AB: only hessian of Y_0
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}


# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters_rescale_identity <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute psi sample
  # parameters$A[,] <- matrix(rnorm(prod(dim(parameters$A))), nrow(parameters$A), ncol(parameters$A))/5
  # browser()
  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )

  # TODO: rescale before or after computing the psi_sim??
  psi_sim <- compute_psi_star_rescaled(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=F, rescale.target=NULL)
  #update covZ
  psi_sam <- compute_psi_star_rescaled(Y, X, Z, parameters, families, Miss, compute_hessian=T, rescale.target=NULL)

  parameters$covZ <- 0.9*parameters$covZ + .1 * cov(psi_sim$Z)*(nrow(psi_sim$Z)/(nrow(psi_sim$Z) - 1))

  # Update AB: only hessian of Y_0
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}




# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters_rescale_biased <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute psi sample
  # parameters$A[,] <- matrix(rnorm(prod(dim(parameters$A))), nrow(parameters$A), ncol(parameters$A))/5
  # browser()
  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )
  # update covZ
  psi_sim <- compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=F)
  parameters$covZ <- 0.9*parameters$covZ + .1 * cov(psi_sim$Z)*(nrow(psi_sim$Z)/(nrow(psi_sim$Z) - 1))

  # # compute Zstar on sample
  # psi_sam <- compute_psi_star(Y, X, Z, parameters, families, Miss, compute_hessian=F)
  # Z <- psi_sam$Z
  #
  # # rescale Z and A
  # resc <- rescale(Z, A, target.cov = parameters$covZ)
  # parameters$A <- resc$A
  # Z <- resc$Z

  # compute psi on this
  psi_sam <- compute_psi_star_rescaled(Y, X, Z, parameters, families, Miss, compute_hessian=T, rescale.target = parameters$covZ)

  # update
  # Update AB: only hessian of Y_0
  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}

# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters_new <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute zhat on sample
  Z0 <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=Z, Miss=Miss)$Zstar
  # rescale
  resc <- rescale(Z0, parameters$A, target.cov=parameters$covZ)
  parameters$A <- resc$A
  Z0 <- resc$Z
  rm(resc)

  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )
  # Obtain Zh
  Zh <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar


  # rescale Zh
  resc <- rescale(Zh, parameters$A, target.cov = parameters$covZ)

  # update covZ
  parameters$covZ <- 0.9*parameters$covZ + .1 * cov(Zh)*(nrow(Zh)/(nrow(Zh) - 1))

  Zh <- resc$Z
  rm(resc)

  # compute psi on this
  psi_sam <- compute_psi_star_known_Z(Y, X, Z0, parameters, families, Miss, compute_hessian=T)
  psi_sim <- compute_psi_star_known_Z(Y_sim$Y, X, Zh, parameters, families, Miss, compute_hessian=F)

  # update

  # Update AB: only hessian of Y_0
  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}


# Update the parameters
# Z is Z.start
#' @parem exponential smoothing coefficient
#' @return list(Z, parameters)
update_parameters_new_2 <- function(Y, X, Z, parameters, families, Miss, alpha, beta) {
  # Compute zhat on sample
  Z0 <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=Z, Miss=Miss)$Zstar
  # rescale
  resc <- rescale(Z0, parameters$A, target.cov=parameters$covZ)
  parameters$A <- resc$A
  Z0 <- resc$Z
  rm(resc)

  # Compute psi simulated
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y),
    Miss = Miss
  )
  # Obtain Zh
  Zh <- compute_zstar(Y_sim$Y, parameters$A, parameters$phi, X, parameters$B, families, start=Y_sim$Z, Miss=Miss)$Zstar


  # rescale Zh
  resc <- rescale(Zh, parameters$A, target.cov = parameters$covZ)

  # update covZ
  parameters$covZ <- 0.9*parameters$covZ + .1 * cov(Zh)*(nrow(Zh)/(nrow(Zh) - 1))

  Zh <- resc$Z
  rm(resc)

  # compute psi on this
  psi_sam <- compute_psi_star(Y, X, Z0, parameters, families, Miss, compute_hessian=T)
  psi_sim <- compute_psi_star(Y_sim$Y, X, Zh, parameters, families, Miss, compute_hessian=F)


  # Update AB: only hessian of Y_0
  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)
  AB_update <- compute_hessian_x_psi(psi_sam$psi_AB - psi_sim$psi_AB, psi_sam$psi_AB_hessian)

  AB_update <- AB_separate(AB_update, ncol(parameters$A))

  phi_update <- psi_sam$psi_phi - psi_sim$psi_phi

  parameters$A <-   parameters$A -   alpha * AB_update$A
  parameters$B <-   parameters$B -   alpha * AB_update$B
  parameters$phi <- parameters$phi - alpha * phi_update #TODO: should be + without the hessian no?

  c(list(Z=psi_sam$Z), parameters)
}




ZX_join <- function(Z, X) {
  if (is.null(X)) {
    Z
  } else {
    cbind(Z,X)
  }
}

AB_separate <- function(AB, q, k=NULL) {
  if (is.null(k)) {
    k <- ncol(AB) - q
  } else {
    stopifnot(q+k == ncol(AB))
  }
  if ( k == 0) {
    A = AB
    B = NULL
  } else {
    A = AB[, 1:q, drop = F]
    B = AB[, (q+1):ncol(AB), drop=F]
  }
  list(A = A, B = B)
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
compute_psi_star <- function(Y, X, Z, parameters, families, Miss, compute_hessian, Z.maxit=100) {
  # overhead computations
  if(!is.null(parameters$B)){
    XB <- X %*% t(parameters$B)
  } else {
    XB <- NULL
  }
  # Imputing step: compute Zstar
  Z <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=Z, Miss=Miss, verbose=F, maxit=Z.maxit)$Zstar

  linpar <- compute_linpar(Z, parameters$A, XB=XB)
  linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)

  ZX <- ZX_join(Z, X)
  #compute psi_star_AB
  psi_AB <- compute_psi_AB(Y, ZX, parameters$phi, linpar_bprime, Miss=Miss)
  psi_phi <- rep(0, length(parameters$phi)) # TODO...

  if(compute_hessian) {
    psi_AB_hessian <- compute_psi_AB_hessian(ZX, parameters$phi, linpar_bprimeprime, Miss)
  } else {
    psi_AB_hessian <- NULL
  }

  list(psi_AB = psi_AB, psi_phi = psi_phi, psi_AB_hessian = psi_AB_hessian, linpar=linpar, Z=Z)
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
  psi_phi <- rep(0, length(parameters$phi)) # TODO...

  if(compute_hessian) {
    psi_AB_hessian <- compute_psi_AB_hessian(ZX, parameters$phi, linpar_bprimeprime, Miss)
  } else {
    psi_AB_hessian <- NULL
  }

  list(psi_AB = psi_AB, psi_phi = psi_phi, psi_AB_hessian = psi_AB_hessian, linpar=linpar, Z=Z)
}


#This returns a list of high-level updates for the parameter
compute_psi_star_rescaled <- function(Y, X, Z, parameters, families, Miss, compute_hessian, Z.maxit=100, rescale.target=NULL) {
  # overhead computations
  if(!is.null(parameters$B)){
    XB <- X %*% t(parameters$B)
  } else {
    XB <- NULL
  }
  # Imputing step: compute Zstar
  Z0 <- compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, start=Z, Miss=Miss, verbose=F, maxit=Z.maxit)$Zstar

  # rescale Zstar to target
  resc <- rescale(Z0, parameters$A, target.cov=rescale.target)
  Z <- resc$Z

  linpar <- compute_linpar(Z, parameters$A, XB=XB)
  linpar_bprime <- compute_linpar_bprime(linpar$linpar, families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar$linpar, families)

  ZX <- ZX_join(Z, X)
  #compute psi_star_AB
  psi_AB <- compute_psi_AB(Y, ZX, parameters$phi, linpar_bprime, Miss=Miss)
  psi_phi <- rep(0, length(parameters$phi)) # TODO...

  if(compute_hessian) {
    psi_AB_hessian <- compute_psi_AB_hessian(ZX, parameters$phi, linpar_bprimeprime, Miss)
  } else {
    psi_AB_hessian <- NULL
  }

  list(psi_AB = psi_AB, psi_phi = psi_phi, psi_AB_hessian = psi_AB_hessian, linpar=linpar, Z=Z0)
}


compute_hessian_x_psi <- function(AB_psi, AB_hessian) {
  prod <- sapply(seq_along(AB_hessian), function(j) {
    as.vector(solve(AB_hessian[[j]], AB_psi[j,]))
  }, simplify=F)
  do.call(rbind, prod)
}

# #not used
# compute_psi_AB_j_hessian <- function(ZX,  phi_j, linpar_bprimeprime_j) {
#   -(t(ZX) %*% (ZX*(linpar_bprimeprime_j)))/phi_j
# }

compute_psi_AB <- function(Y, ZX, phi, linpar_bprime, Miss){
  if (!is.null(Miss)) {
    Y[Miss] <- 0

  }
  (t(Y - linpar_bprime)/phi) %*% ZX
}


compute_psi_AB_hessian <- function (ZX, phi, linpar_bprimeprime, Miss) {
  if (!is.null(Miss)) {
    linpar_bprimeprime[Miss] <- 0
    # compute a list of all hessians
    sapply(1:length(phi), function(j) {
      -(t(ZX) %*% (ZX*(linpar_bprimeprime[,j])))/phi[j]
    }, simplify=F)

  } else {
    # compute a list of all hessians
    sapply(1:length(phi), function(j) {
      -(t(ZX) %*% (ZX*(linpar_bprimeprime[,j])))/phi[j]
    }, simplify=F)

  }
}



# This is only for testing purposes, supposed to be equal to the above
compute_psi_AB_test <- function(Y, ZX, phi, linpar_bprime, Miss){
  if(!is.null(Miss)) {
    Y[miss] <- 0
  }
  diag(1/phi) %*% t(Y - linpar_bprime) %*% ZX
}

compute_psi_AB_j <- function(Y_j, ZX, phi_j, linpar_bprime_j){
  (t(ZX) %*% (Y_j - linpar_bprime_j))/phi_j
}

compute_psi_phi <- function () {

}

compute_AB_initial_values <- function() {
  #consider running the glm below
}

compute_A_glm <- function(Y, Z, X, families, maxit=100) {
  ZX <- ZX_join(Z, X)
  fit <- sapply(1:ncol(Y), function(j) {
    glm(Y[,j] ~0 + ZX, family = families$vec[j], control=list(maxit=maxit))$coef
  }, simplify=F)
  do.call(rbind, fit)
}

compute_error <- function(A1, A2) {
  if(any(dim(A1)!=dim(A2))) stop("Dimensions unequal.")
  A1 <- psych::Procrustes(A1, A2)$loadings
  norm(A1 - A2, type="F")/(norm(A1, type="F") + norm(A2, type="F"))
}

# Update experiments
# ---------


if(0) {
  devtools::load_all()
  set.seed(14251)
  poisson  <- 0
  gaussian <- 0
  binomial <- 20
  q <- 1
  p <- poisson + gaussian + binomial
  fg <- gen_fastgllvm(nobs=300, p=p, q=q, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)), k=1, intercept=T, miss.prob = 0, scale=1)


  zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, Miss=Miss))
  plot(fg$Z, zhat$Zstar)

  #TODO: problem with missing value when q=1

  # zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, linpar$XB, families, Miss=Miss))
  allpar <- initialize_parameters(fg, rescale=F)
  allpar$covZ <- cov(allpar$Z)

  par(mfrow=c(2,1))
  plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=1)
  plot(fg$Z, allpar$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
  par(mfrow=c(1,1))

  # continue with safe stuff

  for(i in 1:50) {
    A.old <- allpar$A
    allpar <- update_parameters_rescale_biased(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=.5/sqrt(i), beta=0)
    print(compute_error(allpar$A, fg$parameters$A))
    if(i%%5 == 1){
    # if(T){
      Sys.sleep(.5)
      par(mfrow=c(2,1))
      plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-5,5), xlim=c(-3,3)); abline(0,1,col=1)
      plot(fg$Z, allpar$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
      par(mfrow=c(1,1))
    }
  }

  # and now finish it!

  for(i in 1:50){
    A.old <- allpar$A
    allpar <- update_parameters_new(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=.5/sqrt(i), beta=0)
    # allpar <- update_parameters_new_2(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=.5/sqrt(i), beta=0)
    # allpar <- update_parameters_rescale(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=0.5/sqrt(i), beta=0)

    # those next are bad but rescaled_biased is super stable
    # allpar <- update_parameters_rescale_biased(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=1/sqrt(i), beta=0)

    error <- norm(A.old-allpar$A, type="F")/(norm(A.old, type="F") + norm(allpar$A, type="F"))

    if(i%%5 == 1){
      Sys.sleep(.5)
      par(mfrow=c(2,1))
      plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-5,5), xlim=c(-3,3)); abline(0,1,col=1)
      plot(fg$Z, allpar$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
      par(mfrow=c(1,1))
    }
    # resc <- rescale(allpar$Z, allpar$A)
    # allpar$A <- resc$A
    # allpar$Z <- resc$Z
    # cat("\nChange:", compute_error(allpar$A, A.old))
    # cat("\nError:", compute_error(allpar$A, fg$parameters$A))
    print(compute_error(allpar$A, fg$parameters$A))
  }
  #
  # library(gllvm)
  # fit.gllvm <- gllvm(y=fg$Y, x=fg$X, formula=~0+x, num.lv=1, family="binomial", sd.errors=F)
  #
  library(gmf)

  # fit.gmf <- gmf(fg$Y, X=fg$X, family=binomial(), p = q)
  # fit.gmf <- gmf(fg$Y, X=fg$X, family=poisson(), p = q)

  library(mirtjml)

  fit.mirtjml <- mirtjml_expr(fg$Y, K= q, tol = 1e-2)
  compute_error(allpar$A, fg$parameters$A)
  compute_error(fit.gmf$v, fg$parameters$A)
  compute_error(fit.mirtjml$A_hat, fg$parameters$A)

  plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-5,5), xlim=c(-3,3))
  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=2); abline(0,1,col=1)
  points(fg$parameters$A, psych::Procrustes(fit.mirtjml$A_hat, fg$parameters$A)$loadings, col=3); abline(0,1,col=1)

  plot(fg$parameters$B, allpar$B)
  points(fg$parameters$B, fit.mirtjml$d_hat, col=3)

}


# TESTS
# --------
if(0) {

  devtools::load_all()
  set.seed(1423)
  poisson  <- 0
  gaussian <- 0
  binomial <- 100
  p <- poisson + gaussian + binomial
  fg <- gen_fastgllvm(nobs=1000, p=p, q=1, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)), k=1, intercept=F, miss.prob = 0)
  zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families, Miss=Miss))
  plot(fg$Z, zhat$Zstar)

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
