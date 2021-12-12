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
#
compute_gradient <- function(){

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
expected_gradient_bernoulli_A <- function(X, par){
  # TODO: this should use a faster version to generate the bernoulli RV
  # TODO: We should not generate the Bernoullis, simply the Z and given Z, obtain the expected gradient directly
  n <- nrow(X)
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
