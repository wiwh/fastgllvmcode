library(tidyverse)
library(psych)
# load all sims
files <- list.files("experiments/GLLVM_binary/low_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/low_dimensions/simres/", file)))


# Simulation settings:
rep <- 50
q <- 2
p.list <- c(20, 40)
n.list <- c(100, 200, 500, 1000)
settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gllvm", "gmf", "mirtjml", "ltm")


# extract all estimated loadings for each simulation for a given setting
# returns a list of loadings matrices, one for each method, and one for the true values
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
compute_errors <- function(par_est, par_true) {
  bias2 <- mean((t(colMeans(par_est) - t(par_true)))^2)
  var   <- mean(colMeans(scale(par_est, scale=F)^2))
  mse <- mean((par_est - par_true)^2)
  mpe <- mean(sqrt(rowSums((par_est-par_true)^2)))/ncol(par_est)
  bias <- colMeans(par_est - par_true)
  list(bias2=bias2, bias=bias, var=var, mse=mse, mpe=mpe)
}

# Bias and variance as function of loading value
ma <- function(series, h=1) {
  n <- length(series)
  ma <- rep(0, n-h)
  for (i in 0:h) {
    ma <- ma + series[(i+1):(n-(h-i))]
  }
  ma/(h+1)
}

trim <- function(matrix, h) {
  apply(matrix, 2, function(series) {
    series[series > h] <- h
    series[series < -h] <- -h
    series
  })
}

# extract the errors
loadings <- sapply(1:nrow(settings), function(i) {
  loadings <- extract_loadings(simres, n=settings[i,1], p = settings[i,2])
  loadings <- lapply(methods, function(method){
    mse <- ma(colMeans((trim(loadings[[method]] - loadings$true, 3)^2)), 5)
    bias <- ma(colMeans(trim(loadings[[method]] - loadings$true, 3))^2, 5)
    mse <- mse-bias

    tibble(mpe=sqrt(mse), bias=bias, true=ma(loadings$true[1,], 5)) %>% mutate(n=settings[i,1], p=settings[i,2], Estimator=method) %>%
      pivot_longer(-c(n, p, Estimator, true), values_to="Values", names_to="Errors")
  })

  do.call(rbind, loadings)
}, simplify=F)
loadings <- do.call(rbind, loadings)


loadings[(loadings$n==1000) & (loadings$Errors=="bias") & (loadings$Estimator=="sprime") & (loadings$p==40), "Values"] <- loadings[(loadings$n==1000) & (loadings$Errors=="bias") & (loadings$Estimator=="sprime") & (loadings$p==40), "Values"]#  + 0.05
loadings <- loadings %>% mutate(Errors = factor(Errors, labels=c("var", "square bias"), levels=c("mpe", "bias")))

# extrat errors due

loadings %>% filter(p==40) %>% filter(Estimator %in% c("prime", "sprime",  "ltm", "mirtjml", "gmf", "gllvm")) %>%
  mutate(n=factor(n, labels=paste("n =", unique(n)), levels=unique(n)))%>%
  ggplot(aes(x=true, y=Values, col=Estimator)) +
  geom_line(linewidth=.8) +
  facet_grid(Errors~n, scale="free_y") +
  theme_bw() +
  xlab("True value of the loadings") +
  ylab("") +
  theme(text=element_text(size=18))


ggsave(file="experiments/GLLVM_binary/low_dimensions/binary_lowdims_bias-variance.png", width=13, height=7)
