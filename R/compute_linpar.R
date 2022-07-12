# Computes the linear parameter of the model given all its elements.
compute_linpar <- function(Z, X, A, B){
  if(dim(X)[2] == 0){
    linpar <- Z %*% t(A)
  } else {
    linpar <- Z %*% t(A) + X %*% t(B)
  }
  linpar
}
