library(tidyverse)
devtools::load_all()
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



rep <- 50
q <- 2
p.list <- c(20, 40)
n.list <- c(100, 200, 500)
settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

i<-1
boxplot(MPE.sims %>% filter(p==settings[i,2] & n == settings[i,1]) %>% dplyr::select(prime, sprime, gmf, ltm), outline=F)
i <- i+1
