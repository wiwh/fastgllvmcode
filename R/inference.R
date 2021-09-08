# The following code is written with readability in mind, not computational efficiency.
# Tremendous computational speedups can be obtained by, for instance, vectorizing the methods.

compute_hessian_estimator <- function(score, theta, Y, c, seed=NULL){
  p <- length(theta)
  perturbation <- sample(c(-c,c), p, replace = T)

  # sometimes, the score requires seeds to compute
  if(!is.null(seed)){
    set.seed(seed)
    cat("\n seed set to ", seed)
  }
  Gk <- score(theta  + perturbation, Y)
  if(!is.null(seed)) set.seed(seed)
  Gk <- Gk - score(theta - perturbation, Y)
  Hk <- Gk %*% t(1/perturbation)
  Hk <- 1/2 * (Hk/2 + t(Hk/2))
  Hk
}
compute_expected_hessian <- function(generator, score, theta, N, c=1e-4, score.seed=F){
  # Returns an estimate of the expectation of the hessian
  # Generator is a function that takes theta and returns a matrix whose rows are the realizations from the true model.
  # Score is a function that takes theta and an observation Y (a vector), and returns a p-vector
  # theta is a p-vector of parameters
  # N is the number of summands to compute the final hessian.
  # score.seed is whether the two scores must be computed with the same seed.

  # Extract information from the arguments
  p <- length(theta)

  # Generate the realizations
  Y <- generator(theta, N)

  H <- matrix(0, p, p)
  # The sum is computed sequentially to reduce the memory requirements
  for(i in 1:N){
    if(score.seed) seed <- 12312+i else seed <- NULL
    H <- H + compute_hessian_estimator(score, theta, Y[i,], c=c, seed=seed)
  }
  H <- H/N
  list(H = H)
}


# here are the functions  used in the document:


# conditional density of Y (n * p matrix) given Z (n * q matrix), returns an n-vector

fYZ2 <- function(Y, Z, par){
  natpar <- Z %*% t(par$A)
  n <- nrow(Y)
  sapply(1:n, function(i){
    prod(exp((Y[i,] * natpar[i, ] - bfunc(natpar[i,]))))
  })
}

dlogfYZ2 <- function(Y, Z, par){
  n <- nrow(Y)
  p <- nrow(par$A)
  q <- ncol(par$A)
  natpar <- Z %*% t(par$A)
  t(sapply(1:n, function(i){
    Zimat <- matrix(Z[i,], nrow=p, ncol=q, byrow=T)
    matderiv <- (Y[i,] - sigmoid(natpar[i,])) * Zimat
    vecderiv <- as.vector(t(matderiv))
  }))
}

