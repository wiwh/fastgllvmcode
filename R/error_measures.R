
#' Compute the mean procrustes error between two loading matrices.
#' @param A1: matrix of estimated loadings
#' @param A2: true matrix of loadings
MPE <- function(A1, A2, rotate=T) {
  if(any(dim(A1)!=dim(A2))) stop("Dimensions unequal.")
  if(rotate) A1 <- psych::Procrustes(A1, A2)$loadings
  norm(A1 - A2, type="F")/prod(dim(A1))
}

#' Compute the mean squared error
#' @param B1: matrix of estimated coefficients
#' @param B2: true matrix of coefficients
MSE <- function(B1, B2) {
  norm(B1 - B2, type="F")/prod(B1)
}

#' Compute the Mean Deviance of a fitted model
#'
#' @param mu: matrix of estimated means
#' @inheritParams fastgllvm

MD <- function(Y, mu, family) {
  families <- generate_families(family, ncol(Y))
  deviance <- Y
  for (family in seq_along(unique(families$vec))) {
    id <- families$id[[family]]
    dev.resids <- families$objects[[family]]$dev.resids
    deviance[,id] <- sapply(families$id[[family]], function(j) {
      dev.resids(Y[,j], mu[,j], 1)
    })
  }
  mean(deviance)
}
