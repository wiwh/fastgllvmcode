library(tidyverse)
library(psych)
# load all sims
files <- list.files("experiments/GLLVM_binary/low_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/low_dimensions/simres/", file)))


# Simulation settings:
rep <- 50
q <- 2
p.list <- c(20, 40)
n.list <- c(50, 100, 500)
settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gllvm", "gmf", "mirtjml", "ltm")


extract_loadings <- function(simres, n, p) {
  simres_subset <- list()

  for (sim in simres) {
    if (sim$model$dimensions$n == n & sim$model$dimensions$p == p) {
      simres_subset <- c(simres_subset, list(sim))
    }
  }

  if(length(simres_subset)==0) return(0)

  o <- order(simres_subset[[1]]$model$parameters$A)
  loadings <- sapply(methods, function(method){
    res <- lapply(simres_subset, function(sim){
      c(as.vector(psych::Procrustes(sim[[method]]$A, sim$model$parameters$A)$loadings))[o]
    })
    do.call(rbind, res)
  }, simplify=F)

  loadings_true <- t(sapply(simres_subset, function(sim){
      c(as.vector(sim$model$parameters$A))[o]
    }))

  c(loadings, list(true = loadings_true))
}



# both par_est and par_true must be matrices
compute_mse <- function(par_est, par_true) {
  bias2 <- mean((t(colMeans(par_est) - t(par_true)))^2)
  var   <- mean(colMeans(scale(par_est, scale=F)^2))
  mse <- mean(colMeans((par_est - par_true)^2))

  c(bias2=bias2, var=var, mse=mse)
}



sapply(1:nrow(settings), function(i) {
  loadings <- extract_loadings(simres, n=settings[i,1], p = settings[i,2])
  if(!is.list(loadings)) return(0)
  sapply(methods, function(method) {
    compute_mse(loadings[[method]], loadings$true)
  })
})



loadings <- extract_loadings(simres, 50, 20)

boxplot(loadings$ltm, outline=F)
boxplot(loadings$sprime, outline=F)
points(1:(ncol(loadings$gllvm)), loadings$true[1,], col=2)
