library(tidyverse)
load("./experiments/high_dimensions/high_dimensions_results.Rdata")
# re-organize the results from setting > iterations > methods to setting > method > iterations

methods  <- c("sa", "sp", "gmf")

i <- 6

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
  q <- sapply(bigsim, function(sim)sim[[1]]$q)
  settings <- cbind(n, p, q)
  settings
}

extract_error <- function(bigsim, sim_num, estimator){
  map(bigsim[[sim_num]], function(mylist) mylist[[estimator]]$error[1]) %>% unlist()
}

# Get the table of median of the mean MSE
mse.median <- t(sapply(1:6, function(i){
  sapply(methods, function(est){
    round(median(extract_error(bigsim, i, est)),2)
  })
}))

mse.mean <- t(sapply(1:6, function(i){
  sapply(methods, function(est){
    round(mean(extract_error(bigsim, i, est)),2)
  })
}))

# combine the two
table.mse <- sapply(seq_along(methods), function(i){
  apply(cbind(mse.mean[,i], mse.median[,i]),1, function(res)paste0(res[1]," (",res[2], ")" ))
})

settings <- extract_settings(bigsim) %>% as.data.frame()

table.mse <- cbind("", settings$p, settings$n, table.mse)

(table.mse.tex <- knitr::kable(table.mse, format="latex", booktabs = T))


# get the table of all timees
table.times <- t(sapply(1:6, function(i){
  round(sapply(methods, function(est) extract_time(bigsim, i, est)),2)
}))

table.times <- cbind("", settings$p, table.times)
(table.times.tex <- knitr::kable(table.times, format="latex", booktabs=T))

