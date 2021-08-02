#' Generates data from a GLLVM model.
#'
#' Nothing to say
#'
#' @param n the number of observations to draw
#' @param mu mean
#' @param sigma sigma
#'
#' @return gllvm a list
#' @export

gen_gllvm <- function(n, mu, sigma){
  stats::rnorm(n, mu, sigma)
}
