compute_phi <- function(fg, return_fastgllvm=F) {
  if (!is.null(fg$families$id$gaussian)) {
    if (!is.null(fg$Miss)) {
      denominator <- fg$dimensions$n - colSums(fg$Miss) - fg$dimensions$q
      warning("compute_phi needs to be improved")
    } else {
      denominator <- fg$dimensions$n - fg$dimensions$q
    }
    fg$parameters$phi[fg$families$id$gaussian] <- colSums((fg$Y[,fg$families$id$gaussian, drop=F] - fg$mean[,fg$families$id$gaussian,drop=F])^2)/denominator
  }
  fg$parameters$phi
}
