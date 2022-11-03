AB_separate <- function(AB, dimensions) {
   A = AB[, 1:dimensions$q, drop = F]
   B = AB[, (dimensions$q+1):ncol(AB), drop=F]
  list(A = A, B = B)
}

get_compute_gradients <- function(method) {
  if (method == "simple") return(compute_gradients_simple)
  if (method == "full") return(compute_gradients_full)
}
