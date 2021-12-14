compute_K <- function(par){
  A <- par$A
  p <- nrow(A)
  q <- ncol(A)
  Psi <- par$Psi
  Ap <- A/Psi
  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
}

#TODO : this must follow the recommendations
compute_learningRate <- function(iter, start, end){
  start*(exp(log(end/start)/iter))**(1:iter) - end
}

# TODO: write this to make it clearer. Compue the gradient specifically for the Binary Data.
# Then we can simply use any method.
bernoulli.gradient <- function(Y, X, par){
  # this should be everything needed to compute the gradient.
}

# Estimate Bernoulli with Sample Paths method
# TODO: add a no-intercept option
bernoulli.estimate.sp <- function(Y, q, X=matrix(1, nrow(Y), 1), H=1, reps=4, iter=250, par.init=NULL, compute.Q=F, verbose=T){
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
      B <- t(Y) %*% X / nobs
      Psi <- rep(1, p)
    }

    A.hist <- matrix(0, iter, p*q)
    B.hist <- matrix(0, iter, p*k)
    crit2.hist <- crit.hist <- rep(0, iter)
    n.Z <- nobs
    dir.cumul <- dir <- rep(0, p*q)

    stepsize.beg <- 10
    stepsize.end <- 9
    rho0 <- stepsize.beg*(exp(log(stepsize.end/stepsize.beg)/iter))**(1:iter) - stepsize.end
    rho0 <- rho0**0 + 10

    Ys <- scale(Y, scale=F)
    Z <- matrix(rnorm(n.Z*q), n.Z, q)
    for(j in 1:reps){
      rho <- rho0 * .7**j
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
          Yprobs <- bernoulli.mean(theta=A, B=B, X=X, Z=Z, seed=h+12312)
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
          A <- A + rho[i] *(SxK %*% Q - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
        } else {
          Q <- diag(q)
          A <- A + rho[i] *(SxK - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
        }
        # now update B
        B <- B + rho[i]/10 * t(Y - Yprobs) %*% X/nobs# %*% XXneg
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
    }
    list(Y=Y, A=A, B=B)
}




# Estimate bernoulli with Stochastic Approximation
bernoulli.estimate.ffa <- function(Y, q, X=matrix(1, nrow(Y), 1), reps=4, iter=250, A.init=NULL, B.init=NULL, Psi=rep(1, ncol(Y)), compute.Q=F, verbose=T){
  # get some parameter values
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

      # Batch gradient
      if(reps==1){
        s <- sample(nobs, max(min(nobs, i/iter * 1000),min(5*q, nobs)), repl=F)
      } else {
        s <- sample(nobs, max(min(nobs, (((1:reps)/reps)**2 * 1000)[j]), min(5*q, nobs)), repl=F)
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
        A <- A + rho[i] *(SxK %*% Q - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
      } else {
        Q <- diag(q)
        A <- A + rho[i] *(SxK - ESxK) # Q should be around Identity in expectation, so it doesn't appear in the second term
      }
      # now update B
      B <- B + rho[i]/10 * t(Y[s,] - Yprobs) %*% X[s,]/n# %*% XXneg
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
  }
  list(Y=Y, A=A, B=B)
}





# Return the mean of the Bernoulli Random Variables.
# Code that efficiently!
bernoulli.mean <- function(theta, X, B, Z=NULL, seed=NULL, n=NULL){
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
