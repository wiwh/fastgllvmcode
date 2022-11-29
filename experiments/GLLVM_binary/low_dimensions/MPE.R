library(tidyverse)
# load all sims
files <- list.files("experiments/GLLVM_binary/low_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/low_dimensions/simres/", file)))

# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gllvm", "gmf", "mirtjml", "ltm")

MPE.sims.list <- lapply(simres, function(sim){
  A0 <- sim$model$parameters$A
  B0 <- sim$model$parameters$B

  res <- sapply(methods, function(method) {
    MPE(sim[[method]]$A, A0)
  })
  names(res) <- methods

  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q)
})

MPE.sims <- do.call(rbind, MPE.sims.list) %>% as_tibble()

boxplot(MPE.sims %>% filter(p==20 & n == 100) %>% dplyr::select(ltm, prime, sprime), outline=F)
