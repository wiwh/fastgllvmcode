#' Computes the linear parameter
#'
#' @inheritParams compute_psi
compute_linpar <- function(Z, A, X, B=NULL, XB=NULL) {
  ZA <- Z %*% t(A)
  if (is.null(B) & is.null(XB)) {
    linpar <- ZA
  } else {
    if(is.null(XB)) XB <- X %*% t(B)
    linpar <- ZA + XB
  }
  list(linpar=linpar, XB=XB, ZA=ZA)
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


compute_mean <- function(fg, linpar=NULL, mean=TRUE) {
  if (is.null(linpar)) {
    fg$linpar <- with(fg, compute_linpar(Z, parameters$A, X, parameters$B))$linpar
  } else {
    fg$linpar <- linpar
  }

  if (mean) {
    fg$mean <- compute_linpar_bprime(fg$linpar, fg$families)
  }
  fg
}
