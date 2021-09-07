# The following code is written with readability in mind, not computational efficiency.
# Tremendous computational speedups can be obtained by, for instance, vectorizing the methods.

compute_hessian_estimator <- function(score, theta, Y, c){
  p <- length(theta)
  perturbation <- sample(c(-c,c), p, replace = T)
  Gk <- score(theta  + perturbation, Y) - score(theta - perturbation, Y)
  Hk <- Gk %*% t(1/perturbation)
  Hk <- 1/2 * (Hk/2 + t(Hk/2))
  Hk
}
compute_expected_hessian <- function(generator, score, theta, N, c=1e-4){
  # Returns an estimate of the expectation of the hessian
  # Generator is a function that takes theta and returns a matrix whose rows are the realizations from the true model.
  # Score is a function that takes theta and an observation Y (a vector), and returns a p-vector
  # theta is a p-vector of parameters
  # N is the number of summands to compute the final hessian.

  # Extract information from the arguments
  p <- length(theta)

  # Generate the realizations
  Y <- generator(theta, N)

  H <- matrix(0, p, p)
  # The sum is computed sequentially to reduce the memory requirements
  for(i in 1:N){
    H <- H + compute_hessian_estimator(score, theta, Y[i,], c=c)
  }
  H <- H/N
  list(H = H)
}
