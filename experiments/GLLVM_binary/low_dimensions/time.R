library(tidyverse)
# load all sims
files <- list.files("experiments/GLLVM_binary/low_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/low_dimensions/simres/", file)))

# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gllvm", "gmf", "mirtjml", "ltm")

time.sims.list <- lapply(simres, function(sim){
  res <- sapply(methods, function(method) {
    sim[[method]]$time[3]
  })
  names(res) <- methods
  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q)
})

time.sims <- do.call(rbind, time.sims.list) %>% as_tibble()

boxplot(time.sims %>% select(all_of(methods)))
