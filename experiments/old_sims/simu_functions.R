
extract_loadings <- function(bigsim, sim_num, estimator){
  loadings <- map(bigsim[[sim_num]], function(mylist) mylist[[estimator]]$A %>% as.vector())
  loadings <- do.call(rbind, loadings)
  loadings
}

extract_time <- function(bigsim, sim_num, estimator){
  map(bigsim[[sim_num]], function(mylist) mylist[[estimator]]$time[1]) %>% unlist() %>% mean()
}

extract_settings <- function(bigsim){
  n <- sapply(bigsim, function(sim)sim[[1]]$n)
  p <- sapply(bigsim, function(sim)sim[[1]]$p)
  settings <- cbind(p, n)
  settings
}
