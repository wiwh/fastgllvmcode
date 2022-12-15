# Compute the mean deviance as a function of n and p for all methods.
devtools::load_all()
library(tidyverse)
# load all sims
files <- list.files("experiments/GLLVM_binary/large_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/large_dimensions/simres/", file)))

# Simulation settings: (same as the simulation file)
q <- 5
p.list <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
n.list <- c(500)
setting <- c("A", "B")


# We compare 5 methods:

settings <- expand.grid(n.list, p.list, setting)
colnames(settings) <- c("n", "p", "setting")
settings$setting <- as.character(settings$setting)


# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gmf", "mirtjml")



# Mean procrustes errors for the loadings
MPE.sims.list <- lapply(simres, function(sim){
  A0 <- sim$model$parameters$A

  res <- sapply(methods, function(method) {
    (MPE(sim[[method]]$A, A0))
  })

  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q, setting=sim$setting)
})

MPE.sims <- do.call(rbind, MPE.sims.list) %>% as_tibble() %>% mutate(across(c(p, all_of(methods)), ~as.numeric(.)))


# Mean squared errors for the intercepts
MSE.sims.list <- lapply(simres, function(sim){
  B0 <- sim$model$parameters$B

  res <- sapply(methods, function(method) {
    MPE(sim[[method]]$B, B0)
  })
  names(res) <- methods

  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q, setting=sim$setting)
})


MSE.sims <- do.call(rbind, MSE.sims.list) %>% as_tibble() %>% mutate(across(c(p, all_of(methods)), ~as.numeric(.)))


# MPE TABLE
# Create the table (we took the median instead of the mean because gllvm is otherwise too bad. For the rest, it does not change the results much.)
table_MPE <- MPE.sims %>% pivot_longer(all_of(methods), values_to="MPE", names_to="Estimator") %>% group_by(setting, p, n, Estimator) %>% summarize(across(MPE,  ~round(mean(.),4))) %>%
  pivot_wider(values_from="MPE", names_from="Estimator") %>% ungroup() %>% dplyr::select(setting, p, all_of(methods)) %>% arrange(setting, p)

knitr::kable(table_MPE %>% filter(setting=="A"), format = "latex")
knitr::kable(table_MPE %>% filter(setting=="B"), format = "latex")
# MSE TABLE
table_MSE <- MSE.sims %>% pivot_longer(all_of(methods), values_to="MPE", names_to="Estimator") %>% group_by(setting, p, n, Estimator) %>% summarize(across(MPE,  ~round(mean(.),4))) %>%
  pivot_wider(values_from="MPE", names_from="Estimator") %>% ungroup() %>% dplyr::select(setting, p, all_of(methods)) %>% arrange(setting, p)

knitr::kable(table_MSE %>% filter(setting=="A"), format = "latex")
knitr::kable(table_MSE %>% filter(setting=="B"), format = "latex")

