compute_error <- function(A1, A2, rotate=T) {
  if(any(dim(A1)!=dim(A2))) stop("Dimensions unequal.")
  if(rotate) A1 <- psych::Procrustes(A1, A2)$loadings
  norm(A1 - A2, type="F")/(norm(A1, type="F") + norm(A2, type="F"))
}
