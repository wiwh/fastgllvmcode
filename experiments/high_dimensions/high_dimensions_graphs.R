library(tidyverse)
load("./experiments/estimators_comparison/estimators_comparison_results.Rdata")
# re-organize the results from setting > iterations > methods to setting > method > iterations

methods  <- c("sa", "sp","gllvm", "gmf")


i <- 3

gllvm <- map(bigsim[[i]], function(mylist) mylist$gllvm$A %>% as.vector())
gllvm <- do.call(rbind, gllvm)
gllvm.time <- map(bigsim[[i]], function(mylist) mylist$gllvm$time[1]) %>% unlist() %>% mean()

sa <- map(bigsim[[i]], function(mylist) mylist$sa$A %>% as.vector())
sa <- do.call(rbind, sa)
sa.time <- map(bigsim[[i]], function(mylist) mylist$sa$time[1]) %>% unlist() %>% mean()

sp <- map(bigsim[[i]], function(mylist) mylist$sp$A %>% as.vector())
sp <- do.call(rbind, sp)
sp.time <- map(bigsim[[i]], function(mylist) mylist$sp$time[1]) %>% unlist() %>% mean()

gmf <- map(bigsim[[i]], function(mylist) mylist$gmf$A %>% as.vector())
gmf <- do.call(rbind, gmf)
gmf.time <- map(bigsim[[i]], function(mylist) mylist$gmf$time[1]) %>% unlist() %>% mean()

p <- bigsim[[i]][[i]]$p
true_loadings <- c(seq(2,-2, l=p), seq(-1, 1, l=p))

boxplot(sa, outline=F, main=sa.time)
points(true_loadings, col=2)

boxplot(sp, outline=F, main=sp.time)
points(true_loadings, col=2)

boxplot(gmf, outline=F, main=gmf.time)
points(true_loadings, col=2)

boxplot(gllvm, outline=F, main=gllvm.time)
points(true_loadings, col=2)
