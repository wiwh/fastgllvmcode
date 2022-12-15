# Compute the mean deviance as a function of n and p for all methods.
devtools::load_all()
library(tidyverse)
# load all sims
files <- list.files("experiments/GLLVM_poisson/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_poisson/simres/", file)))

# Simulation settings:
q <- 5
p.list <- c(100, 200, 200, 300, 400, 500)
n.list <- c(100, 500)
setting <- c("A", "B")



# We compare 5 methods:

settings <- expand.grid(n.list, p.list, setting)
colnames(settings) <- c("n", "p", "setting")
settings$setting <- as.character(settings$setting)

settings <- settings[1:5,]


# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gmf")


# extract all estimated loadings for each simulation for a given setting
# returns a list of loadings matrices, one for each method, and one for the true values
extract_loadings <- function(simres, n, p, setting) {
  simres_subset <- list()

  for (sim in simres) {
    if (!is.list(sim)) next()
    if (sim$model$dimensions$n == n & sim$model$dimensions$p == p & sim$setting == setting) {
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
  loadings <- extract_loadings(simres, n=settings[i,1], p = settings[i,2], setting=settings[i,3])
  loadings <- lapply(methods, function(method){
    mse <- colMeans((trim(loadings[[method]] - loadings$true, 3)^2))
    bias <- colMeans(loadings[[method]] - loadings$true)

    tibble(mpe=sqrt(mse), bias=bias, true=loadings$true[1,]) %>% mutate(n=settings[i,1], p=settings[i,2], setting=settings[i,3], Estimator=method) %>%
      pivot_longer(-c(n, p, setting, Estimator, true), values_to="Values", names_to="Errors")
  })

  do.call(rbind, loadings)
}, simplify=F)
loadings <- do.call(rbind, loadings)


# loadings[(loadings$n==1000) & (loadings$Errors=="bias") & (loadings$Estimator=="sprime") & (loadings$p==40), "Values"] <- loadings[(loadings$n==1000) & (loadings$Errors=="bias") & (loadings$Estimator=="sprime") & (loadings$p==40), "Values"] - .05
loadings <- loadings %>% mutate(Errors = factor(Errors, labels=c("MSE", "bias"), levels=c("mpe", "bias")),
                                setting = factor(setting, labels=paste("Setting", c("A", "B")), levels=c("A", "B")))

# extrat errors due

pA <- loadings %>% filter(setting=="Setting A", p %in% c(100, 200)) %>%
  mutate(p = factor(p, labels=paste("p =", (1:10)*100), levels=(1:10)*100)) %>%
  ggplot(aes(x=true, y=Values, col=Estimator)) +
  geom_smooth(size=.8, se=F) +
  facet_grid(Errors~p, scale="free_y") +
  theme_bw() +
  xlab("True value of the loadings") +
  ylab("") +
  ggtitle("Setting A") +
  theme(text=element_text(size=18))

pA


pB <- loadings %>% filter(setting=="Setting B", p %in% c(100, 500, 1000)) %>%
  mutate(p = factor(p, labels=paste("p =", (1:10)*100), levels=(1:10)*100)) %>%
  ggplot(aes(x=true, y=Values, col=Estimator)) +
  geom_smooth(size=.8, se=F) +
  facet_grid(Errors~p, scale="free_y") +
  theme_bw() +
  xlab("True value of the loadings") +
  ylab("") +
  ggtitle("Setting B") +
  theme(text=element_text(size=18))

library(ggpubr)
ggarrange(pA, pB, nrow=1, common.legend = TRUE, legend="bottom")

ggsave(file="experiments/GLLVM_binary/large_dimensions/binary_largedims_settingAB_bias-variance.png", width=13, height=7)
