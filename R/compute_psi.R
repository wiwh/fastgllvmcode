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

  # Update AB: each his own hessian
  # psi_sim <- compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=T)
  # AB_update <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian) -
  #              compute_hessian_x_psi(psi_sim$psi_AB, psi_sim$psi_AB_hessian)
  # AB_update <- AB_separate(AB_update, ncol(parameters$A))

  # Update AB: only hessian of Y_0
  psi_sim <- compute_psi_star(Y_sim$Y, X, Y_sim$Z, parameters, families, Miss, compute_hessian=F)
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
compute_psi_star <- function(Y, X, Z, parameters, families, Miss, compute_hessian) {
  # overhead computations
  if(!is.null(parameters$B)){
    XB <- X %*% t(parameters$B)
  } else {
    XB <- NULL
  }
  # Imputing step: compute Zstar
  Z <- compute_zstar(Y, parameters$A, parameters$phi, XB, families, start=Z, Miss=Miss, verbose=F)$Zstar

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
  A1 <- psych::Procrustes(A1, A2)$loadings
  norm(A1 - A2, type="F")/(norm(A1, type="F") + norm(A2, type="F"))
}

# Update experiments
# ---------


if(0) {
  devtools::load_all()
  set.seed(1423)
  fg <- gen_fastgllvm(nobs=1000, p=200, q=1, family=c(rep("poisson", 200), rep("gaussian", 0), rep("binomial", 00)), k=5, intercept=F, miss.prob = 0)
  # zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, linpar$XB, families, Miss=Miss))
  allpar <- allpar.init <- initialize_parameters(fg)
  # Zstar <- compute_zstar(fg$Y, allpar$A, allpar$phi, fg$X %*% t(allpar$B), fg$families, start=allpar$Z)$Zstar

  # par(mfrow=c(2,1))
  # plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=1)
  # plot(fg$Z, allpar$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
  # par(mfrow=c(1,1))
  # allpar <- c(list(Z=fg$Z), fg$parameters)


  for(i in 1:1000){
    A.old <- allpar$A
    allpar <- update_parameters(fg$Y, fg$X, allpar$Z, allpar, fg$families, fg$Miss, alpha=1, beta=0)

    error <- norm(A.old-allpar$A, type="F")/(norm(A.old, type="F") + norm(allpar$A, type="F"))
    if(1){
      if(i%%5 == 1){
      # if(T){
        Sys.sleep(1)
        par(mfrow=c(2,1))
        plot(fg$parameters$A, psych::Procrustes(allpar$A, fg$parameters$A)$loadings, ylim=c(-5,5), xlim=c(-3,3)); abline(0,1,col=1)
        plot(fg$Z, allpar$Z, ylim=c(-3,3), xlim=c(-3,3)); abline(0,1,col=2); abline(0,-1,col=2)
        par(mfrow=c(1,1))
      }
    }
    # resc <- rescale(allpar$Z, allpar$A)
    # allpar$A <- resc$A
    # allpar$Z <- resc$Z
    cat("\nChange:", compute_error(allpar$A, A.old))
    cat("\nError:", compute_error(allpar$A, fg$parameters$A))
  }
}

# TESTS
# --------
if(0) {

  devtools::load_all()
  set.seed(1423)
  fg <- gen_fastgllvm(nobs=200, p=50, q=1, family=c(rep("poisson", 00), rep("gaussian", 000), rep("binomial", 50)), k=1, intercept=F, miss.prob = 0)
  zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, linpar$XB, families, Miss=Miss))
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
