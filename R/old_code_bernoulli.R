# This is some old code to be re-factored.

gllvm.rbinary <- function(n, A, Z=NULL, B=matrix(0, nrow(A), 1), X=matrix(1, n, 1)){
  p <- nrow(A)
  q <- ncol(A)
  # if Z is not provided, generate it
  if(is.null(Z)) Z <- matrix(rnorm(n*q), n, q)

  linpar <- Z %*% t(A) + X %*% t(B)
  # if the uniform sample is not provided, generate it
  # TODO: test that this is strictly equivalent with the same seed (it is, but test it)
  # Warning! this is extremely inefficient. Provide sunif if possible.
  Y <- sapply(1:p, function(j){
    rbinom(n = n,
           size = 1,
           prob = 1/(1 + exp(-linpar[,j])))
  })
  list(Y=Y, Z=Z, A=A, B=B, X=X)
}


# X is the history of X
generator.bernoulli <- function(n, theta, Z=NULL, sunif=NULL, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  p <- nrow(theta)
  q <- ncol(theta)
  # if Z is not provided, generate it
  if(is.null(Z)) Z <- matrix(rnorm(n*q), n, q)

  linpar <- Z %*% t(theta)
  # if the uniform sample is not provided, generate it
  # TODO: test that this is strictly equivalent with the same seed (it is, but test it)
  # Warning! this is extremely inefficient. Provide sunif if possible.
  if(is.null(sunif)){
    X <- sapply(1:p, function(j){
      rbinom(n = n,
             size = 1,
             prob = 1/(1 + exp(-linpar[,j])))
    })
    # else use the uniform sample
  } else {
    X <- 1* (linpar > sunif)
  }
  list(X=X, Z=Z)
}

# compute the expectation of the covariance given z
bernoulli.cov.exp <- function(theta, Z=NULL, seed=NULL, n=NULL){
  if(is.null(Z) && is.null(n)) stop("Provide a value for Z or n.")
  q <- ncol(theta)
  p <- nrow(theta)
  if(exists(".Random.seed")){
    INITIALSEED <- .Random.seed
    lockBinding("INITIALSEED", environment())
    on.exit(.Random.seed <<- INITIALSEED)
  }
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Z)){
    Z <- matrix(rnorm(n*q), n, q)
  } else {
    n <- nrow(Z)
  }
  linpar <- Z %*% t(theta)
  # matrices of conditional expectations (probabilities)
  probs <- 1/(1+exp(-linpar))
  probs.zeromean <- scale(probs, scale=F)
  M <- t(probs.zeromean) %*% probs.zeromean/n
  Mdiag <- colMeans(probs); Mdiag <- Mdiag*(1-Mdiag)
  diag(M) <- Mdiag
  M
}

bernoulli.cov.exp.covariates <- function(theta, X, B, Z=NULL, seed=NULL, n=NULL){
  if(is.null(Z) && is.null(n)) stop("Provide a value for Z or n.")
  q <- ncol(theta)
  p <- nrow(theta)
  if(exists(".Random.seed")){
    INITIALSEED <- .Random.seed
    lockBinding("INITIALSEED", environment())
    on.exit(.Random.seed <<- INITIALSEED)
  }
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Z)){
    Z <- matrix(rnorm(n*q), n, q)
  } else {
    n <- nrow(Z)
  }
  linpar <- Z %*% t(theta) + X %*% t(B)
  # matrices of conditional expectations (probabilities)
  probs <- 1/(1+exp(-linpar))
  probs.zeromean <- scale(probs, scale=F)
  M <- t(probs.zeromean) %*% probs.zeromean/n
  Mdiag <- colMeans(probs); Mdiag <- Mdiag*(1-Mdiag)
  diag(M) <- Mdiag
  M
}

bernoulli.probs <- function(theta, Z=NULL, seed=NULL, n=NULL){
  if(is.null(Z) && is.null(n)) stop("Provide a value for Z or n.")
  q <- ncol(theta)
  p <- nrow(theta)
  if(exists(".Random.seed")){
    INITIALSEED <- .Random.seed
    lockBinding("INITIALSEED", environment())
    on.exit(.Random.seed <<- INITIALSEED)
  }
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Z)){
    Z <- matrix(rnorm(n*q), n, q)
  } else {
    n <- nrow(Z)
  }
  linpar <- Z %*% t(theta)
  # matrices of conditional expectations (probabilities)
  1/(1+exp(-linpar))
}
bernoulli.probs.covariates <- function(theta, X, B, Z=NULL, seed=NULL, n=NULL){
  if(is.null(Z) && is.null(n)) stop("Provide a value for Z or n.")
  q <- ncol(theta)
  p <- nrow(theta)
  if(exists(".Random.seed")){
    INITIALSEED <- .Random.seed
    lockBinding("INITIALSEED", environment())
    on.exit(.Random.seed <<- INITIALSEED)
  }
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Z)){
    Z <- matrix(rnorm(n*q), n, q)
  } else {
    n <- nrow(Z)
  }
  linpar <- Z %*% t(theta) + X %*% t(B)
  # matrices of conditional expectations (probabilities)
  1/(1+exp(-linpar))
}

bernoulli.ffa.stochastic <- function(Y, q, Psi=rep(1, ncol(Y))){
  nobs <- nrow(Y)
  p <- ncol(Y)
  iter <- 200
  reps <- 5
  A <- diag(1, p, q)
  A.hist <- matrix(0, iter, p*q)
  crit2.hist <- crit.hist <- rep(0, iter)
  n.Z <- 5*nobs
  Z <- matrix(rnorm(n.Z*q), n.Z, q)
  dir.cumul <- dir <- rep(0, p*q)

  stepsize.beg <- 400
  stepsize.end <- 10
  rho0 <- stepsize.beg*(exp(log(stepsize.end/stepsize.beg)/iter))**(1:iter) - stepsize.end

  Ys <- scale(Y, scale=F)
  for(j in 1:reps){
    rho <- rho0 * .7**j
    A.old <- A
    for(i in 1:iter){
      # Corrected gradients
      Ap <- A/Psi
      K  <-  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)

      s <- sample(nobs, max(min(nobs, (((1:reps)/reps)**2 * 500)[j]),q), repl=F)
      # This is for the ffa step
      SxK <- t(Ys[s,]) %*% (Ys[s,] %*% t(K/(length(s)-1)))  # maybe divide earlier by sqrt(nrow(Y))?

      # And this is for the expectation of the ffa step
      Yp <- bernoulli.probs(A, Z[s,])
      Yp <- scale(Yp, scale=F)
      diagYY <- attr(Yp, "scaled:center"); diagYY <- diagYY * (1-diagYY)
      vars <- colMeans(Yp^2)
      ESxK <- t(Yp) %*% (Yp %*% t(K/(length(s)-1)))
      # HERE I NEED TO MODIFy THE EFFECT OF HAVING DIAGONAL ELEMENTS OF YY not EQUAL TO THOSE OF T(Y)Y
      ESxK <- ESxK + t(K) * (diagYY - vars)
      # SxK is now the expected covariance times t(K). Time to compute the step. We need Q
      Q <- backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K %*% ESxK)
      A <- A + rho[i] *(SxK %*% Q - ESxK) # %*% Q # Q should be around Identity in expectation, so it doesn't appear in the second term
      A.hist[i, ] <- as.vector(A)
      # crit.hist[i] <- norm(ESxK - SxK)
      crit2.hist[i] <- norm(Q-diag(q))
      cat("\ni: ", i, " - norm:", crit2.hist[i])
    }
  }
  list(Y=Y, A=A)
}


bernoulli.ffa.covariate <- function(Y, q, X=matrix(1, nrow(Y), 1), reps=4, iter=250, A.init=NULL, B.init=NULL, Psi=rep(1, ncol(Y)), compute.Q=T, verbose=T){
  nobs <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)
  if(!is.null(A.init)) A <- A.init else A <- diag(1, p, q)
  if(!is.null(B.init)) B <- B.init else B <- t(Y) %*% X / nobs
  A.hist <- matrix(0, iter, p*q)
  B.hist <- matrix(0, iter, p*k)
  crit2.hist <- crit.hist <- rep(0, iter)
  n.Z <- nobs
  dir.cumul <- dir <- rep(0, p*q)

  stepsize.beg <- 200
  stepsize.end <- 1
  rho0 <- stepsize.beg*(exp(log(stepsize.end/stepsize.beg)/iter))**(1:iter) - stepsize.end

  Ys <- scale(Y, scale=F)
  Z <- matrix(rnorm(n.Z*q), n.Z, q)
  for(j in 1:reps){
    rho <- rho0 * .7**j
    A.old <- A
    for(i in 1:iter){
      # Corrected gradients
      Ap <- A/Psi
      K  <-  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)

      if(reps==1){
        s <- sample(nobs, max(min(nobs, i/iter * 1000),min(5*q, nobs)), repl=F)
      } else {
        s <- sample(nobs, max(min(nobs, (((1:reps)/reps)**2 * 1000)[j]), min(5*q, nobs)), repl=F)
      }

      # This is for the ffa step
      SxK <- t(Ys[s,]) %*% (Ys[s,] %*% t(K/(length(s)-1)))  # maybe divide earlier by sqrt(nrow(Y))?

      # And this is for the expectation of the ffa step
      Yprobs <- bernoulli.probs.covariates(theta=A, B=B, X=X[s,], Z=Z[s,, drop=F])
      Yp <- scale(Yprobs, scale=F)
      diagYY <- attr(Yp, "scaled:center"); diagYY <- diagYY * (1-diagYY)
      vars <- colMeans(Yp^2)
      ESxK <- t(Yp) %*% (Yp %*% t(K/(length(s)-1)))
      # HERE I NEED TO MODIFy THE EFFECT OF HAVING DIAGONAL ELEMENTS OF YY not EQUAL TO THOSE OF T(Y)Y
      ESxK <- ESxK + t(K) * (diagYY - vars)
      # SxK is now the expected covariance times t(K). Time to compute the step. We need Q
      if(compute.Q){
        Q <- backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K %*% ESxK)
      } else {
        Q <- diag(1, q)
      }
      A <- A + rho[i] *(SxK %*% Q - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
      # now update B
      B <- B + rho[i]/10 * t(Y[s,] - Yprobs) %*% X[s,]/n# %*% XXneg
      # save
      B.hist[i, ] <- as.vector(B)
      A.hist[i, ] <- as.vector(A)
      # crit.hist[i] <- norm(ESxK - SxK)
      crit2.hist[i] <- norm(Q-diag(q))
      if(verbose)cat("\ni: ", i, " - norm:", crit2.hist[i])
    }
  }
  list(Y=Y, A=A, B=B)
}


generator <- function(t, seed=1243){
  rbinom(nobs,1,sigmoid(rnorm(nobs, 0,1),c(0,t)))
}
estimator <- function(x, seed=2141, maxiter=10){
  for(i in 1:maxiter){
    set.seed(seed)
    # given ah, estimate y, linear model uesh yo
    yh <- ah[2]/(ah[2]+1)*(x-mean(x))
    # modify it so variance match expectation
    sdx <- get.sdx(ah)
    yh <- yh/sd(yh) * sdx
    # do normal one
    (ah <- c(mean(x), mean((x-mean(x))*yh/sdx**2)))
  }
  ah
}

# tgrid <- (1:100)/10
# sims <- sapply(tgrid, function(t) estimator(generator(t))[2])
# plot(tgrid, sims)


# MULTI DIMENSIONS
# ---------------
ffa.cov <- function(covX, q, maxiter=100, eps=1e-4, A.init = NULL, savepath=F, Psi.fix=NULL){
  p <- ncol(covX)
  if(savepath)path <- list()
  if(!is.null(Psi.fix)){
    if(eps != 0) warning("Psi.fix was not NULL; eps=0 was used, maxiter is used as unique ending criterion.")
    eps <- 0
  }
  if(!is.null(A.init)) A <- A.init else A <- diag(1, p,q)
  if(!is.null(Psi.fix)) Psi <- Psi.fix else Psi <- rep(1, p)
  for(i in 1:maxiter){
    Psi.old <- Psi
    Ap <- A/Psi # storing this yields to HUGE IMPROVEMENTS
    K <- solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
    SxK <- covX %*% t(K)
    A <- SxK %*% (backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K%*%A))
    if(!is.null(Psi.fix)){
      Psi <- Psi.fix
    } else {
      Psi <- diag(covX - A %*% t(A))
    }
    Psi[Psi < 1e-10] <- 1e-10
    if(savepath) path[[i]] <-  list(A=A, Psi=Psi)
    if(eps && (mean(abs(Psi.old-Psi)))<eps) break
  }
  list(L=A, com=Psi, niter=i, path=if(savepath)path)
}

ffa <- function(X, q, maxiter=100, eps=1e-4, A.init = NULL, Psi.fix = NULL, savepath=F, compute.Q=T){
  if(savepath)path <- list()
  n <- nrow(X)
  p <- ncol(X)
  if(!is.null(Psi.fix)){
    if(eps != 0) warning("Psi.fix was not NULL; eps=0 was used, maxiter is used as unique ending criterion.")
    eps <- 0
  }
  if(!is.null(Psi.fix)) Psi <- Psi.fix else Psi <- rep(1, p)
  if(!is.null(A.init)) A <- A.init else A <- diag(1, p,q)
  X <- scale(X, scale=F)
  vars <- colMeans(X^2)
  for(i in 1:maxiter){
    Psi.old <- Psi
    Ap <- A/Psi # storing this yields to HUGE IMPROVEMENTS
    K <- solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
    SxK <- t(X) %*% (X %*% t(K/n))
    if(compute.Q){
      A <- SxK %*% (backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K%*%A))

    } else {
      A <- SxK
    }
    if(!is.null(Psi.fix)){
      Psi <- Psi.fix
    } else{
      Psi <- vars - as.vector(sapply(1:p, function(j) sum(A[j,]^2)))

    }
    Psi[Psi < 1e-10] <- 1e-10
    if(savepath) path[[i]] <-  list(A=A, Psi=Psi)
    if(eps && (mean(abs(Psi.old-Psi)))<eps) break
  }
  list(L=A, com=Psi, niter=i, path=if(savepath)path)
}

# X  = probs is logit(linpar), used to compute the covariance matrix.
# the offdiagonal elements correspond to those of cov(probs)
# however, the diagaonal elements do not, so there is extra code to change the effect of this
ffa.probs <- function(X, q, maxiter=100, eps=1e-4, A.init = NULL, Psi.fix = NULL, savepath=F, compute.Q = T){
  if(savepath)path <- list()
  p <- ncol(X)
  n <- nrow(X)

  if(!is.null(Psi.fix)){
    if(eps != 0) warning("Psi.fix was not NULL; eps=0 was used, maxiter is used as unique ending criterion.")
    eps <- 0
  }
  if(!is.null(A.init)) A <- A.init else A <- diag(1, p,q)
  if(!is.null(Psi.fix)) Psi <- Psi.fix else Psi <- rep(1, p)

  X <- scale(X, scale=F)

  diagXX <- attr(X, "scaled:center"); diagXX <- diagXX * (1-diagXX)
  vars <- colMeans(X^2)

  for(i in 1:maxiter){
    Psi.old <- Psi
    Ap <- A/Psi # storing this yields to HUGE IMPROVEMENTS
    K <- solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
    SxK <- t(X) %*% (X %*% t(K/n))
    # HERE I NEED TO MODIFY THE EFFECT OF HAVING DIAGONAL ELEMENTS OF XX not EQUAL TO THOSE OF T(X)X

    SxK <- SxK + t(K) * (diagXX - vars)

    if(compute.Q){
      A <- SxK %*% (backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K%*%A))

    } else {
      A <- SxK
    }

    if(!is.null(Psi.fix)){
      Psi <- Psi.fix
    } else {
      Psi <- diagXX - as.vector(sapply(1:p, function(j) sum(A[j,]^2)))
    }
    Psi[Psi < 1e-10] <- 1e-10
    if(savepath) path[[i]] <-  list(A=A, Psi=Psi)
    if(eps && (mean(abs(Psi.old-Psi)))<eps) break
  }
  list(L=A, com=Psi, niter=i, path=if(savepath)path)
}
