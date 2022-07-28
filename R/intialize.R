library(fastfactoranalysis)

initialize_parameters <- function(fastgllvm, target=NULL, rescale=T) {
  cat("\nInitializing...\n")
  pb <- txtProgressBar(style=3, label="Initializing", width=40)
  with(fastgllvm, {
    # TODO: try initializing using gllvm
    Y.transformed <- Y
    if(length(families$id$poisson) > 0) {
      Y.transformed[,families$id$poisson] <- log(Y.transformed[,families$id$poisson] + 0.1) * 1.5
    }

    if(length(families$id$binomial) > 0) {
      Y.transformed[,families$id$binomial] <- (Y.transformed[,families$id$binomial] - 0.5) *6
    }

    NAs <- is.na(Y)
    Y.transformed[NAs] <- 0

    setTxtProgressBar(pb, 0.2,  title="Initializing B")
    # Initialize XB
    if(dimensions$k >0) {
      B <- t(lm(Y.transformed ~ 0+X)$coef) * (1+colMeans(NAs)) # TODO: ffa should initialize with missing values as well
      XB <- X %*% t(B)
      Y.transformed <- Y.transformed - XB
    } else {
      B <- NULL
      XB <- NULL
    }

    setTxtProgressBar(pb, 0.4, title="Initializing A")
    # Initialize A and phi
    fit.ffa <- ffa(Y.transformed, dimensions$q, maxiter=10) # TODO: there must be some Z here as well...
    A <- fit.ffa$loadings * (1+colMeans(NAs))  # TODO: make ffa work with NAs, but this trick is ok for initialization
    phi <- rep(1, dimensions$p)
    phi[families$id$gaussian] <- fit.ffa$communalites[families$id$gaussian]

    # go to target
    if (!is.null(target)) {
      A <- psych::Procrustes(A, target)$loadings
    }
    setTxtProgressBar(pb, 0.9, title="Initializing Z")
    # initialize Z
    Z <- init_Z(Y.transformed, A, phi)
    if(rescale) {
      rescaled <- rescale(Z, A)
      A=rescaled$A
      Z=rescaled$Z
    }
    setTxtProgressBar(pb, 1, title="Initialization Complete.")
    close(pb)
    # re-scale zhat and A
    list(A=A, B=B, phi=phi, Z=Z)
  })
}

init_Z <- function(Y, A, phi) {
  Y %*% (A/phi) %*% solve(t(A) %*% (A/phi))
}
# rescale Z to have unit diagonal variance, and A so that ZA remains the same value
rescale <- function(Z, A) {
  #TODO: rescale B as well here? if Z is not centered, recenter it and change B so that XB + ZA remains the same value
  Z <- scale(Z, scale=F)
  C <- chol((t(Z) %*% Z)/nrow(Z))
  Z <- Z %*% solve(C)
  A <- A %*% t(C)

  # TODO: TEST THAT: Zt(A) (after centing Z) remains the same
  list(Z=Z, A=A)
}

if(0) {
  devtools::load_all()
  set.seed(1231)
  fg <- gen_fastgllvm(nobs=1000, p=1000, q=20, family=c(rep("poisson", 500), rep("gaussian", 0), rep("binomial", 500)), k=1, intercept=1, miss.prob = 0.5)
  fg <- gen_fastgllvm(nobs=1000, p=100, q=1, family=c(rep("poisson", 0), rep("gaussian", 0), rep("binomial", 100)), k=1, intercept=F, miss.prob = 0)
  init <- initialize_parameters(fg, target=fg$parameters$A, rescale=F)

  plot(fg$parameters$A, init$A); abline(0,1,col=2)
  points(fg$Z, init$Z, col=2)
  points(fg$parameters$phi, init$phi, col=3)
  points(fg$parameters$B, init$B, col=4)
}
