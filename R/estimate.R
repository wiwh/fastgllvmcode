get_K <- function(A, B, phi){
  p <- nrow(A)
  q <- ncol(A)
  Ap <- A/phi
  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
}

#TODO : this must follow the recommendations
get_learning_rate.seq <- function(learning_rate.start, learning_rate.end, iter, type="exp", reps=1, reps.decreased_rate=0.5, alpha=0.7){
  if(type=="exp"){
    learning_rate.seq <- rep(0, iter*reps)
    for(j in 1:reps){
      rho <- exp((log(learning_rate.end) - log(learning_rate.start))/(iter-1))
      learning_rate.seq[((j-1)*iter+1):(j*iter)] <- learning_rate.start * rho**(0:(iter-1)) * reps.decreased_rate**(j-1)
    }
    learning_rate.seq
  } else if (type=="sa"){
    learning_rate.seq <- rep(0, iter*reps)
    for(j in 1:reps){
      # alpha <- log(learning_rate.end/learning_rate.start)/log((1+5)/(iter+5))
      rho <- learning_rate.start*(1+5)**alpha/(1:iter + 5)**alpha
      rho <- rho * learning_rate.start / rho[1]
      learning_rate.seq[((j-1)*iter+1):(j*iter)] <- rho * reps.decreased_rate**(j-1)
    }
  }
  learning_rate.seq
}


# Estimate Bernoulli with Sample Paths method
# TODO: add a no-intercept option
bernoulli.estimate.sp <- function(Y, q, X=matrix(1, nrow(Y), 1), H=1, reps=1, iter=250, par.init=NULL, compute.Q=F, verbose=T, learning_rate.start=10, learning_rate.end=1, learning_rate.type="exp", tol=1e-5){
    # get some dimensions
    nobs <- nrow(Y)
    p <- ncol(Y)
    k <- ncol(X)

    # get some parameter values
    if(!is.null(par.init)){
      A <- par.init$A
      B <- par.init$B
      Psi <- par.init$Psi
    } else {
      A <- diag(1, p, q)
      # B <- t(Y) %*% X / nobs
      B <- matrix(0, p, k)
      Psi <- rep(1, p)
    }

    A.hist <- matrix(0, iter, p*q)
    B.hist <- matrix(0, iter, p*k)
    crit2.hist <- crit.hist <- rep(0, iter)
    n.Z <- nobs
    dir.cumul <- dir <- rep(0, p*q)

    Ys <- scale(Y, scale=F)
    Z.list <- lapply(1:H, function(h){
      set.seed(h+12312)
      matrix(rnorm(nobs*q), nobs, q)
    })
    learning_rate.seq <- get_learning_rate.seq(learning_rate.start, learning_rate.end, iter, type=learning_rate.type)
    A.old <- A
    for(i in 1:iter){
      # Corrected gradients
      Ap <- A/Psi
      K  <-  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)

      # This is for the sp step
      SxK <- t(Ys) %*% (Ys %*% t(K/(nobs-1)))  # maybe divide earlier by sqrt(nrow(Y))?

      # And this is for the expectation of the ffa step
      ES.list <- vector(mode="list", length=H)
      for(h in 1:H){
        # always the same seed
        Yprobs <- bernoulli.mean(theta=A, B=B, X=X, Z=Z.list[[h]])
        # Yprobs <- bernoulli.mean(theta=A, B=B, X=X, n=nobs, seed=h+12312)
        Yp <- scale(Yprobs, scale=F)

        diagYY <- attr(Yp, "scaled:center"); diagYY <- diagYY * (1-diagYY)
        vars <- colMeans(Yp^2)
        ESxK <- t(Yp) %*% (Yp %*% t(K/(nobs-1)))
        # HERE I NEED TO MODIFy THE EFFECT OF HAVING DIAGONAL ELEMENTS OF YY not EQUAL TO THOSE OF T(Y)Y
        ESxK <- ESxK + t(K) * (diagYY - vars)
        ES.list[[h]] <- ESxK
      }

      ESxK <- Reduce("+", ES.list)/length(ES.list)

      # ESxK is now the expected covariance times t(K). Time to compute the step. We need Q?
      if(compute.Q){
        Q <- backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K %*% ESxK)
        A <- A + learning_rate.seq[i] *(SxK %*% Q - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
      } else {
        Q <- diag(q)
        A <- A + learning_rate.seq[i] *(SxK - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
      }

      # stop values of loadings that are diverging due to Heywood
      broken.A <- abs(A)>10
      A[broken.A] <- 10 * sign(A[broken.A])
      # now update B
      B <- B + learning_rate.seq[i]/5 * t(Y - Yprobs) %*% X/nobs# %*% XXneg
      # save
      B.hist[i, ] <- as.vector(B)
      A.hist[i, ] <- as.vector(A)
      if(compute.Q){
        crit.hist[i] <- norm(Q-diag(q))
      } else {
        crit.hist[i] <- norm(ESxK - SxK)
      }
      if(verbose)cat("\ni: ", i, " - norm:", crit.hist[i])
      if(crit.hist[i] < tol) break()
    }
    par <- list(A=A,
                B=B,
                ifelse(exists("Psi"), Psi, rep(1, nrow(A))),
                family="bernoulli",
                p=nrow(A),
                q=ncol(A))
    list(Y=Y, A=A, B=B, A.hist=A.hist, B.hist=B.hist, par=par)
}



# batch: False, or a proportion of observations to use every iterations
# Estimate bernoulli with Stochastic Approximation
bernoulli.estimate.ffa <- function(Y, q, X=matrix(1, nrow(Y), 1), iter=250, batch=F, reps=4, reps.decreased_rate=0.7, learning_rate.start=20, learning_rate.end=.1, learning_rate.type="exp", A.init=NULL, B.init=NULL, Psi=rep(1, ncol(Y)), compute.Q=F, verbose=T, tol=1e-5){
  # get some parameter values
  nobs <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(X)
  if(!is.null(A.init)) A <- A.init else A <- diag(1, p, q)
  if(!is.null(B.init)) B <- B.init else B <- matrix(0, p, k) # t(Y) %*% X / nobs

  A.hist <- matrix(0, reps*iter, p*q)
  B.hist <- matrix(0, reps*iter, p*k)
  crit.hist <- rep(0, iter)
  n.Z <- nobs
  dir.cumul <- dir <- rep(0, p*q)


  learning_rate.seq <- get_learning_rate.seq(learning_rate.start, learning_rate.end, iter, reps=reps, reps.decreased_rate = reps.decreased_rate, type=learning_rate.type)

  # standardize Y and initialize parameters
  Ys <- scale(Y, scale=F)

  i <- 1
  while(i < (iter*reps)){
    Z <- matrix(rnorm(n.Z*q), n.Z, q)
    i <- i+1
    # Corrected gradients
    Ap <- A/Psi
    K  <-  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)

    if(!batch){
      s <- 1:nobs
    } else {
      s <- sample(nobs, size=batch*nobs, replace=F)
    }

    # This is for the ffa step
    SxK <- t(Ys[s,]) %*% (Ys[s,] %*% t(K/(length(s)-1)))  # maybe divide earlier by sqrt(nrow(Y))?

    # And this is for the expectation of the ffa step
    Yprobs <- bernoulli.mean(theta=A, B=B, X=X[s,], Z=Z[s,, drop=F])
    Yp <- scale(Yprobs, scale=F)

    diagYY <- attr(Yp, "scaled:center"); diagYY <- diagYY * (1-diagYY)
    vars <- colMeans(Yp^2)
    ESxK <- t(Yp) %*% (Yp %*% t(K/(length(s)-1)))
    # HERE I NEED TO MODIFy THE EFFECT OF HAVING DIAGONAL ELEMENTS OF YY not EQUAL TO THOSE OF T(Y)Y
    ESxK <- ESxK + t(K) * (diagYY - vars)
    # SxK is now the expected covariance times t(K). Time to compute the step. We need Q?
    if(compute.Q){
      Q <- backsolve(chol(K %*%  SxK), diag(1,q)) %*% chol(K %*% ESxK)
      A <- A + learning_rate.seq[i] *(SxK %*% Q - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
    } else {
      Q <- diag(q)
      A <- A + learning_rate.seq[i] *(SxK - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
    }

    # stop values of loadings that are diverging due to Heywood
    broken.A <- abs(A)>10
    A[broken.A] <- 10 * sign(A[broken.A])

    # now update B
    B <- B + learning_rate.seq[i]/5 * t(Y[s,] - Yprobs) %*% X[s,]/nobs# %*% XXneg
    # save
    A.hist[i, ] <- as.vector(A)
    B.hist[i, ] <- as.vector(B)

    if(compute.Q){
      crit.hist[i] <- norm(Q-diag(q))
    } else {
      crit.hist[i] <- (mean((A.hist[i, ] - A.hist[i-1,])^2)/var(A.hist[i,]) + mean((B.hist[i, ] - B.hist[i-1,])^2))/var(B.hist[i,]) # /learning_rate.seq[i]
    }
    if(verbose)cat("\ni: ", i, " - norm:", crit.hist[i], " - learning rate:", learning_rate.seq[i])
    # check if the criterion is small enough to jump to the next "repetition", where the learning rate increases again
    if(i<(iter*reps) && crit.hist[i] < tol){
      # jump to next rep
      reps.jumps <- iter*(1:reps)
      inext <- reps.jumps[which(i < reps.jumps)[1]]
      cat("\n\nnext!!!", inext)

      # fill in the histories
      A.hist[(i+1):inext,] <- matrix(rep(A.hist[i,], inext-i), nrow=(inext-i), byrow = T)
      B.hist[(i+1):inext,] <- matrix(rep(B.hist[i,], inext-i), nrow=(inext-i), byrow = T)
      i <- inext
    }
  }
  par <- list(A=A,
              B=B,
              ifelse(exists("Psi"), Psi, rep(1, nrow(A))),
              family="bernoulli",
              p=nrow(A),
              q=ncol(A),
              k=k)
  list(Y=Y, A=A, B=B, A.hist=A.hist, B.hist=B.hist, crit.hist=crit.hist, par=par)
}

#' Returns the (Monte-Carlo) approximations of $E[Y]$ and $E[Y^\top Y]K^\top$.
#' @param A: the matrix of loadings
#' @param B: the matrix of fixed effect coefficients
#' @param X: the matrix of covariates
#' @param Z: a list of H generated Z samples
#' @param K: the matrix so that $\hat Z = K Y$
#' @param family: an object of class "family"
get_expectations <- function(A, B, X, Z, K, family){
  EY <- lapply(Z, function(Zh)family$linkinv(compute_natpar(A, B, Zh, X)))
  EYYK <- lapply(EY, function(EYh){
    EYh.c <- scale(EYh, scale=F)
    EYYh.diag <- colSums(EYh.c^2) / (nrow(EYh.c))
    EYh.mean <- attr(EYh.c, "scaled:center")
    EYh.var <- family$variance(EYh.mean) * phi  # TODO: check that we need to multiply by phi... or how to model overdispersion here
    EYYKh <- t(EYh.c) %*% (EYh.c %*% t(K/(nrow(EYh.c)-1)))
    # The diagonal of EYZi is wrong. We need to replace its current value (EYYi.diag) by its true value (EYi.var).
    # However, we do not want to compute EYY, so we remove its effect after the multiplication by t(K).
    # t(K) is pre-multiplied because resp.var and EYYi.diag are vectors; it would be post-multiplied if they were
    # diagonal matrices.
    EYYKh <- EYZh + t(K) * (EYh.var - EYYh.diag)
    EYYKh
  })
  list(
    EY = Reduce("+", EY)/length(EY),
    EYYK = Reduce("+", EYZ)/length(EYZ)
  )
}

#' Returns the psi function approximated by H monte carlo samples
#' @param Y: uncentered data
#' @param Y.c: Y that are CENTERED, they must have an attribute named `scaled:center` that represent the removed mean
#' @param A: the matrix of loadings
#' @param B: the matrix of fixed effect coefficients
#' @param phi: the vector of response-specific scale parameters
#' @param X: the matrix of covariates
#' @param family: an object of class "family"
#' @param generate_Z a function returned by generate_Z_functionfactory that returns a list
get_Psi <- function(Y, Y.c, A, B, phi, X, family, generate_Z){
  stopifnot(!is.null(attr(Y.c, "scaled:center")))
  K  <- get_K(A, B, phi)
  Z  <- generate_Z()
  Exp <- get_expectations(A, B, X, Z, K, family)

  Psi <- list(
    A = t(Y.c) %*% (Y.c %*% t(K/(nrow(Y.c)-1))) - Exp$EYYK,
    B = t(Y - Exp$EY) %*% X/nrow(x),
    phi = rep(0, p)  # TODO: implement
  )

  Psi
}

#' Returns a fitted object of class "fastgllvm".
#' @param Y: a matrix of dimensions $n\times q$
#' @param q: the number of latent variables
#' @param X: either 0 (no covariates, no intercept), 1 (no covariates but an intercept), or a matrix of covariates. If the latter and an intercept is desired, it must be included in X as a column of 1.
#' @param family: either one of ("gaussian", "poisson", "binomial"), the function name or the corresponding object of class "family" obtained by calling the function.
#' @param family: one of "SA" or "SP"
#' @param H: how many samples of Z to draw: if $\in (0,1)$, then a batch method is used and only that proportion of data is used to estimate the model.
#'
#' @export
fastgllvm <- function(Y, q=1, X=1, family=binomial(), method="SA", H=1, A.init=NULL, B.init=NULL, phi.init=NULL, iter=250, reps=4, reps.decreased_rate=0.7, learning_rate.start=20, learning_rate.end=.1, learning_rate.type="exp", verbose=T, tol=1e-5){
  stopifnot(is.matrix(Y))
  n <- nrow(Y)
  p <- ncol(Y)

  if(lengnth(X) == 1){
    stopifnot(X == 0 | X == 1)
    X <- matrix(X, n, 1)
  }

  k <- ncol(X)

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }


  if(!is.null(A.init)) A <- A.init else A <- diag(1, p, q)
  if(!is.null(B.init)) B <- B.init else B <- matrix(0, p, k)

  A.hist <- matrix(0, reps*iter, p*q)
  B.hist <- matrix(0, reps*iter, p*k)
  phi.hist <- matrix(0, reps*iter, p)
  crit.hist <- rep(0, iter)
  # dir.cumul <- dir <- rep(0, p*q)


  learning_rate.seq <- get_learning_rate.seq(learning_rate.start, learning_rate.end, iter, reps=reps, reps.decreased_rate = reps.decreased_rate, type=learning_rate.type)

  Y.c <- scale(Y, scale=F)
  generate_Z <- generate_Z_functionfactory(n, q, method=method, H=H)

  i <- 0
  while(i < (iter*reps)){
    i <- i+1
    Psi <- get_Psi(Y, Y.c, A, B, phi, X, family, generate_Z)

    # udate A
    A <- A + learning_rate.seq[i] * Psi$A
    broken.A <- abs(A)>10
    A[broken.A] <- 10 * sign(A[broken.A])

    # update B
    B  <- B + learning_rate.seq[i]/5 * Psi$B

    # update phi
    phi <- phi + learning_rate.seq[i] * Psi$phi

    # save
    A.hist[i, ] <- as.vector(A)
    B.hist[i, ] <- as.vector(B)
    phi.hist[i, ] <- as.vector(B)

    if(verbose)cat("\ni: ", i, " - norm:", crit.hist[i], " - learning rate:", learning_rate.seq[i])
    # check if the criterion is small enough to jump to the next "repetition", where the learning rate increases again
    if(i<(iter*reps) && crit.hist[i] < tol){
      # jump to next rep
      reps.jumps <- iter*(1:reps)
      inext <- reps.jumps[which(i < reps.jumps)[1]]
      cat("\n\nnext!!!", inext)

      # fill in the histories
      A.hist[(i+1):inext,] <- matrix(rep(A.hist[i,], inext-i), nrow=(inext-i), byrow = T)
      B.hist[(i+1):inext,] <- matrix(rep(B.hist[i,], inext-i), nrow=(inext-i), byrow = T)
      i <- inext
    }
  }
  par <- list(A=A,
              B=B,
              ifelse(exists("Psi"), Psi, rep(1, nrow(A))),
              family="bernoulli",
              p=nrow(A),
              q=ncol(A),
              k=k)
  list(Y=Y, A=A, B=B, A.hist=A.hist, B.hist=B.hist, crit.hist=crit.hist, par=par)
}

