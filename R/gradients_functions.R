
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
