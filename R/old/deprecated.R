# deprecated functions

#' Generates parameters for a gllvm
#'
#' returns a list of parameters.
#'
#'
gen_par <- function(p, q, k=0, family="normal", A=NULL, B=NULL, Psi=NULL, scale=0.5){
  warning("gen_par is deprecated. we don't use par anymore.")
  # generate the (p, q) matrix of loadings
  if(is.null(A)) A <- matrix(rnorm(p*q), p, q)*scale
  # generate the vector of scale parameters
  if(is.null(Psi)) Psi <- runif(p, 0.5, 2)
  Psi[family=="bernoulli"] <- 1
  # generate the (p, k) matrix of coefficients
  if(is.null(B)){
    if(k==0){
      B <- matrix(0, p, 1)
    } else {
      B <- matrix(runif(p*k, -1, 1), p, k)
    }
  }
  list(A=A, B=B, Psi=Psi, family=family, p=p, q=q, k=k)
}

#' Modifies par with new parameters
#'
change_par <- function(par, A=NULL, B=NULL, Psi=NULL){
  warning("change_par is deprecated.")
  if(!is.null(A)) par$A <- A
  if(!is.null(B)) par$B <- B
  if(!is.null(Psi)) par$Psi <- Psi
  par
}

sigmoid <- function(x){
  warning("sigmoid is deprecated. use invlink of the respective family.")
  1/(1+exp(-x))
}
