AB_separate <- function(AB, dimensions) {
   A = AB[, 1:dimensions$q, drop = F]
   if (dimensions$k > 0) {
    B = AB[, (dimensions$q+1):ncol(AB), drop=F]
   } else {
     B <- NULL
   }
  list(A = A, B = B)
}

par_to_vec <- function(par_list) {
  AB <- t(cbind(par_list$A, par_list$B))
  c(as.vector(AB), par_list$phi, as.vector(par_list$covZ))
}

vec_to_par <- function(par_vec, dim) {
  AB.last <- dim$p * (dim$q + dim$k)
  phi.last <- AB.last + dim$p
  covZ.last <- length(par_vec)

  AB <- matrix(par_vec[1:AB.last], dim$p, dim$q + dim$k, byrow = T)
  if(dim$k >0 ) {
    A <- AB[,1:dim$q, drop=F]
    B <- AB[,(dim$q+1):ncol(AB), drop=F]
  } else {
    A <- AB
    B <- NULL
  }
  phi <- par_vec[(AB.last+1):phi.last]
  covZ <- matrix(par_vec[(phi.last+1):covZ.last], dim$q, dim$q)
  list(A=A, B=B, phi=phi, covZ=covZ)
}

