# boxplot

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

# extract the quantiles
loadings <- sapply(1:nrow(settings), function(i) {
  loadings <- extract_loadings(simres, n=settings[i,1], p = settings[i,2])
  loadings <- lapply(methods, function(method){
    as_tibble(t(apply(loadings[[method]]- loadings$true, 2, function(x) quantile(x, c(.05, .5, .95))))) %>% mutate(loading=1:n(), true=loadings$true[1,], n=settings[i,1], p=settings[i,2], Estimator=method)
      # pivot_longer(-c(n, p, Estimator, loading), values_to="Quantile", names_to="Prob")
  })

  do.call(rbind, loadings)
}, simplify=F)

loadings <- do.call(rbind, loadings)
loadings %>% filter(n %in% c(500, 1000), p%in%c(40)) %>% filter(Estimator %in% c("prime", "ltm", "gmf", "mirtjml")) %>%
  ggplot(aes(x=true, y=`50%`, col=Estimator)) +
  geom_point() +
  geom_line()+
  geom_ribbon(aes(ymin=`5%`, ymax=`95%`), linetype=2, alpha=0.1) +
  facet_grid(n~p)

ma <- function(series, h=1) {
  n <- length(series)
  ma <- rep(0, n-h)
  for (i in 0:h) {
  ma <- ma + series[(i+1):(n-(h-i))]
  }
  ma/(h+1)
}
# extract the loadings
loadings <- sapply(1:nrow(settings), function(i) {
  loadings <- extract_loadings(simres, n=settings[i,1], p = settings[i,2])
  loadings <- lapply(methods, function(method){
    as_tibble(loadings[[method]]) %>% mutate(loading=1:n(), n=settings[i,1], p=settings[i,2], Estimator=method) %>%
      pivot_longer(-c(n, p, Estimator, loading), values_to="Estimate", names_to="Loading") %>% mutate(true=rep(loadings$true[1,], n()/length(loadings$true[1,])))
  })

  do.call(rbind, loadings)
}, simplify=F)
loadings <- do.call(rbind, loadings)

loadings %>% filter(n ==1000, p==40) %>% filter(Estimator %in% c("prime", "ltm", "gmf", "mirtjml")) %>%
  ggplot(aes(x=true, y=Estimate, col=Estimator)) +
  geom_smooth()+
  geom_abline(col="red")

