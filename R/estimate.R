compute_K <- function(par){
  A <- par$A
  p <- nrow(A)
  q <- ncol(A)
  Psi <- par$Psi
  Ap <- A/Psi
  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
}

# TODO: write this to make it clearer. Compue the gradient specifically for the Binary Data.
# Then we can simply use any method.
bernoulli.gradient <- function(Y, X, par){
  # this should be everything needed to compute the gradient.
}

#TODO : this must follow the recommendations
get_learning_rate.seq <- function(learning_rate.start, learning_rate.end, iter, type="exp", reps=1, reps.decreased_rate=0.7){
  learning_rate.seq <- rep(0, iter*reps)
  for(j in 1:reps){
    rho <- exp((log(learning_rate.end) - log(learning_rate.start))/(iter-1))
    learning_rate.seq[((j-1)*iter+1):(j*iter)] <- learning_rate.start * rho**(0:(iter-1)) * reps.decreased_rate**(j-1)
  }
  learning_rate.seq
}

# Estimate Bernoulli with Sample Paths method
# TODO: add a no-intercept option
bernoulli.estimate.sp <- function(Y, q, X=matrix(1, nrow(Y), 1), H=1, reps=1, iter=250, par.init=NULL, compute.Q=F, verbose=T, learning_rate.start=10, learning_rate.end=1){
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
    learning_rate.seq <- get_learning_rate.seq(learning_rate.start, learning_rate.end, iter)
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
    }
    list(Y=Y, A=A, B=B, A.hist=A.hist, B.hist=B.hist)
}



# batch: False, or a proportion of observations to use every iterations
# Estimate bernoulli with Stochastic Approximation
bernoulli.estimate.ffa <- function(Y, q, X=matrix(1, nrow(Y), 1), iter=250, batch=F, reps=4, reps.decreased_rate=0.7, learning_rate.start=20, learning_rate.end=.1, A.init=NULL, B.init=NULL, Psi=rep(1, ncol(Y)), compute.Q=F, verbose=T, tol=1e-5){
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


  learning_rate.seq <- get_learning_rate.seq(learning_rate.start, learning_rate.end, iter, reps=reps, reps.decreased_rate = reps.decreased_rate)

  # standardize Y and initialize parameters
  Ys <- scale(Y, scale=F)
  Z <- matrix(rnorm(n.Z*q), n.Z, q)

  i <- 1
  while(i < (iter*reps)){
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
  list(Y=Y, A=A, B=B, A.hist=A.hist, B.hist=B.hist, crit.hist=crit.hist)
}





# Return the mean of the Bernoulli Random Variables.
# Code that efficiently!
bernoulli.mean <- function(theta, X, B, Z=NULL, seed=NULL, n=NULL){
  if(is.null(Z) && is.null(n)) stop("Provide a value for Z or n.")
  q <- ncol(theta)
  p <- nrow(theta)
  # if(exists(".Random.seed")){
  #   INITIALSEED <- .Random.seed
  #   lockBinding("INITIALSEED", environment())
  #   on.exit(.Random.seed <<- INITIALSEED)
  # }
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Z)){
    Z <- matrix(rnorm(n*q), n, q)
  } else {
    n <- nrow(Z)
  }
  linpar <- Z %*% t(theta) + X %*% t(B)
  # matrices of conditional expectations (probabilities)
  sigmoid(linpar)
}





gradient_bernoulli_A <- function(Y, X, par){
  n <- nrow(Y)
  K <- compute_K(par)
  Zstar <- Y %*% t(K)
  natpar <- compute_natpar(par, Zstar, X)
  t(Y - sigmoid(natpar)) %*% (Zstar/n)
}

# compute bernoulli probabilities
compute_bernoulli_probabilities <- function(X, par){

  natpar <- compute_natpar(par, Zstar, X)
}


# Estimate the Expected Gradient for the Loadings
expected_gradient_bernoulli_A <- function(X, par, n=nrow(X)){
  # TODO: this should use a faster version to generate the bernoulli RV
  # TODO: We should not generate the Bernoullis, simply the Z and given Z, obtain the expected gradient directly
  Y <- gen_gllvm(n, family="bernoulli", X=X, par=par)$Y
  gradient_bernoulli_A(Y, X, par)
}

corrected_gradient_bernoulli_A <- function(Y, X, par){
  grad <- gradient_bernoulli_A(Y, X, par)
  Egrad <- expected_gradient_bernoulli_A(X, par)
  grad - Egrad
}


profile_score_gllvm <- function(){

}

score_gllvm <- function(params){

}
