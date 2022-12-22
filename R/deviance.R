compute_deviance <- function(fg) {
  deviance <- fg$Y
  for (family in seq_along(unique(fg$families$vec))) {
    id <- fg$families$id[[family]]
    dev.resids <- fg$families$objects[[family]]$dev.resids
    deviance[,id] <- sapply(fg$families$id[[family]], function(j) {
      dev.resids(fg$Y[,j], fg$mean[,j], 1)
    })
  }
  deviance
}
