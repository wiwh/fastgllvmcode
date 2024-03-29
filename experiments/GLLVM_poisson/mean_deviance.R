# Compute the mean deviance as a function of n and p for all methods.
devtools::load_all()
library(tidyverse)
# load all sims
files <- list.files("experiments/GLLVM_poisson/simres/")
simres <- lapply(files, function(file) readRDS(paste0("experiments/GLLVM_poisson/simres/", file)))

# Simulation settings:
q <- 5
p.list <- c(100, 200, 300, 400, 500, 1000)
n.list <- c(100, 500)
setting <- c("A", "B")



# We compare 5 methods:

settings <- expand.grid(n.list, p.list, setting)
colnames(settings) <- c("n", "p", "setting")
settings$setting <- as.character(settings$setting)


# n x p x q x (model ,  5 methods)

methods <- c("prime", "sprime", "gmf")



# compute the mean deviance given a particular setting of n and p
extract_MD <- function(simres, n, p, setting) {
  simres_subset <- list()

  for (sim in simres) {
    if (sim$n == n & sim$p == p & sim$setting == setting) {
      simres_subset <- c(simres_subset, list(sim))
    }
  }

  if(length(simres_subset)==0) return(0)

  sapply(methods, function(method){
    res <- sapply(simres_subset, function(sim){
      # compute mean:
      fg <- sim$model
      fg$parameters$A <- sim[[method]]$A
      fg$parameters$B <- sim[[method]]$B
      fg$parameters$Z <- sim[[method]]$Z
      fg$mean <- compute_mean(fg, return_object = F)$mean
      median(compute_deviance(fg))
    })
  })
}


trim <- function(matrix, h) {
  apply(matrix, 2, function(series) {
    series[series > h] <- h
    series[series < -h] <- -h
    series
  })
}



# extract the errors (this takes a while... so better load the compiled error file)
errors <- sapply(1:nrow(settings), function(i) {
  errors <- extract_MD(simres, n=settings[i,1], p = settings[i,2], setting=settings[i,3])
  errors <- as_tibble(errors) %>% mutate(n=settings[i,1], p=settings[i,2], setting=settings[i,3])

}, simplify=F)
errors <- do.call(rbind, errors)

saveRDS(errors, file="experiments/GLLVM_poisson/mean_deviance.rds")

errors <- readRDS("experiments/GLLVM_poisson/mean_deviance.rds")

e2 <- errors %>% filter(p %in% c(100, 500, 1000)) %>% pivot_longer(-c(n, p, setting), names_to="Estimator", values_to = "MD") %>% mutate(n= paste0("n = ", n), p = factor(p, levels = c(100, 500, 1000), labels= paste0("p = ", c(100, 500, 1000))), setting = paste0("Setting ", setting))

ggplot(e2, aes(x=Estimator, y=MD, color=Estimator)) +
  geom_boxplot() +
  facet_grid(p~setting) +
  ylab("Mean deviance") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none")

ggsave(file="experiments/GLLVM_poisson/poisson_largedims_mean_deviance.png", width=13, height=7)


#
# spline_to_list <- function(series) {
#   sp <- smooth.spline(series, spar=.1)
#   tibble(p=seq(100, 500, l=length(sp$y)), error=sp$y)
# }
#
# sp <- errors %>% pivot_longer(-c(n, p, setting), names_to="Estimator", values_to = "MD") %>% group_by(setting, Estimator) %>% summarize(across(MD, ~spline_to_list(.))) %>% ungroup() %>%
#   mutate(p = MD$p, error = MD$error) %>% mutate(MD=NULL) %>% as_tibble() %>% mutate(setting = factor(setting, levels=c("A", "B"), labels = paste("Setting", c("A", "B"))))
# ggplot(sp, aes(x=p, y=error, col=Estimator)) +
#   facet_grid(.~setting) +
#   geom_line(size=1) +
#   theme_bw() +
#   theme(text = element_text(size=20))+
#   scale_x_continuous()+
#   ylab("Mean Deviance")


# ggsave(file="experiments/GLLVM_binary/large_dimensions/binary_largedims_mean_deviance.png", width=13, height=7)
