#' Computes the linear parameter
#'
#' @inheritParams compute_psi
compute_linpar <- function(Z, A, X, B=NULL, XB=NULL) {
  if (is.null(B) & is.null(XB)) {
    ZA <- Z %*% t(A)
    linpar <- ZA
  } else {
    if(is.null(XB)) XB <- X %*% t(B)
    ZA <- Z %*% t(A)
    linpar <- Z %*% t(A) + XB
  }
  linpar_list <- list(linpar=linpar, XB=XB, ZA=ZA)
  linpar_list
}

compute_linpar_bprime <- function(linpar, families) {
  if(is.matrix(linpar)) {
    for(i in seq_along(families$id)){
      linpar[,families$id[[i]]] <- families$objects[[i]]$linkinv(linpar[,families$id[[i]]])
    }
  } else if (is.vector(linpar)) {
    for(i in seq_along(families$id)){
      linpar[families$id[[i]]] <- families$objects[[i]]$linkinv(linpar[families$id[[i]]])
    }
  } else {
    stop("TypeError linpar")
  }
  linpar
}

compute_linpar_bprimeprime <- function(linpar, families) {
  stopifnot(is.matrix(linpar))
  for(i in seq_along(families$id)){
    linpar[,families$id[[i]]] <- families$objects[[i]]$mu.eta(linpar[,families$id[[i]]])
  }
  linpar
}
