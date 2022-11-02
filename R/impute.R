#' Initialize Imputation
#'
#' @param

initialize_imputation <- function(Y) {
  means <- colMeans(Y, na.rm=T)
  sd  <- apply(Y, 2, sd, na.rm=T)
}
