library(tidyverse)
devtools::load_all()
# load all sims
files <- list.files("experiments/GLLVM_binary/low_dimensions/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_binary/low_dimensions/simres/", file)))

# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gllvm", "gmf", "mirtjml", "ltm")

# Mean procrustes errors for the loadings
MPE.sims.list <- lapply(simres, function(sim){
  A0 <- sim$model$parameters$A

  res <- sapply(methods, function(method) {
    MPE(sim[[method]]$A, A0)
  })
  names(res) <- methods

  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q)
})

MPE.sims <- do.call(rbind, MPE.sims.list) %>% as_tibble()


# Mean squared errors for the intercepts
MSE.sims.list <- lapply(simres, function(sim){
  B0 <- sim$model$parameters$B

  res <- sapply(methods, function(method) {
    MPE(sim[[method]]$B, B0)
  })
  names(res) <- methods

  c(res, n=sim$model$dimensions$n, p=sim$model$dimensions$p, q=sim$model$dimensions$q)
})


MSE.sims <- do.call(rbind, MSE.sims.list) %>% as_tibble()


# MPE TABLE
# Create the table (we took the median instead of the mean because gllvm is otherwise too bad. For the rest, it does not change the results much.)
table_MPE <- MPE.sims %>% pivot_longer(all_of(methods), values_to="MPE", names_to="Estimator") %>% group_by(p, n, Estimator) %>% summarize(across(MPE,  ~round(mean(.),3))) %>%
  pivot_wider(values_from="MPE", names_from="Estimator") %>% dplyr::select(p, n, prime, sprime, ltm, gmf, mirtjml, gllvm)

knitr::kable(table_MPE, format = "latex")
# MSE TABLE
table_MSE <- MSE.sims %>% pivot_longer(all_of(methods), values_to="MPE", names_to="Estimator") %>% group_by(p, n, Estimator) %>% summarize(across(MPE,  ~round(mean(.),3))) %>%
  pivot_wider(values_from="MPE", names_from="Estimator") %>% dplyr::select(p, n, prime, sprime, ltm, gmf, mirtjml, gllvm)

knitr::kable(table_MSE, format = "latex")

