# Constructor
# -----------

#' Generates a fastgllvm object
#'
#' @param A the matrix of loadings.
#' @param B the matrix of fixed effect coefficients, of dimensions p * k
#' @param phi a vector of scale parameters.
#' @param X either 0 (no covariates, no intercept), 1 (an intercept), or a matrix of n * k covariates (with, possibly, the first column of 1s being an intercept)
#'
#'
#' @return a list corresponding to the model
new_fastgllvm <- function(A, B, phi, family, X){
  stopifnot(is.matrix(A))
  stopifnot(is.matrix(B))
  stopifnot(is.vector(phi))
  stopifnot(attr(family, "class")=="family")

  p <- nrow(A)
  q <- ncol(A)
  if(q >= p)
    stop("The number of latent variables (q) must be strictly smaller than the number of observed variables (p).")
  if(nrow(B) != p)
    stop("The number of rows of B must be equal to the number of observed variables (p).")
  if(!is.vector(phi) || length(phi) != p)
    stop("phi must be a vector of length p.")

  if(!is.matrix(X)){
    if(X==0 & length(B) > 0)
      stop("No intercept given, yet fixed effect coefficients are given. Should there be an intercept?")
    if(X==1 & ncol(B) !=1)
      stop("The dimensions of the fixed effect coefficients matrix B is inconsistent with the value of X implying only an intercept is required.")
  } else {
    if(nrow(X) != p || ncol(X) != k)
      stop("X must be  either 0 or 1, or a matrix of dimensions p * k.")
  }

  fastgllvm <- list(A=A,
       B=B,
       phi=phi,
       family=family,
       X=X)
  class(fastgllvm) <- "family"
  fastgllvm
}

# Methods
# -------
predict.fastgllvm <- function(fastgllvm){}
update.fastgllvm <- function(fastgllvm){}
simualte.fastgllvm <- function(fastgllvm){}
plot.fastgllvm <- function(fastgllvm){}
print.fastgllvm <- function(fastgllvm){}
