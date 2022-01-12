library(tidyverse)
library(scales)
library(gridExtra)
library(zoo)
library(latex2exp)

# load results for gllvm
load("./experiments/portion_data/portion_data_results_small_gllvm.Rdata")
bigsim_small_gllvm <- bigsim
rm(bigsim)

# load results for SA, SP, GMF
load("./experiments/portion_data/portion_data_results_small.Rdata")
bigsim_small <- bigsim1
rm(bigsim1)

load("./experiments/portion_data/portion_data_results_large.Rdata")
bigsim_large <- bigsim
rm(bigsim)

# combine small and large sims for SA, SP, GMF
bigsim <- c(bigsim_small, bigsim_large)
rm(bigsim_small, bigsim_large)


# extract the times and error and construct a tibble
types <- c("sa", "sp", "gmf")
times <- map_dfc(types, function(type) map(bigsim, function(sim)sim[[type]]$time[1]) %>% unlist() %>% as.vector)
colnames(times) <- types
times <- map_dfc(times, rollmean, 3)

#number of parameters
numpar <- map(bigsim, function(sim) sim$dat.p*sim$dat.q) %>% unlist() %>% as.vector
error <- map_dfc(types, function(type) map(bigsim, function(sim)sim[[type]]$error[1]) %>% unlist() %>% as.vector) %>%
  mutate_all(function(x)x/(numpar))
colnames(error) <- types
error <- map_dfc(error, rollmean, 3)

# do the same for gllvm
types.gllvm <- c("gllvm")
times.gllvm <- map_dfc(types.gllvm, function(type) map(bigsim_small_gllvm, function(sim)sim[[type]]$time[1]) %>% unlist() %>% as.vector)
colnames(times.gllvm) <- types.gllvm
times.gllvm <- map_dfc(times.gllvm, rollmean, 3)

numpar <- map(bigsim_small_gllvm, function(sim) sim$dat.p*sim$dat.q) %>% unlist() %>% as.vector
error.gllvm <- map_dfc(types.gllvm, function(type) map(bigsim_small_gllvm, function(sim)sim[[type]]$error[1]) %>% unlist() %>% as.vector) %>%
  mutate_all(function(x)x/(numpar))
colnames(error.gllvm) <- types.gllvm
error.gllvm <- map_dfc(error.gllvm, rollmean, 3)

# complete gllvm with na
nrow.tot <- nrow(times)
nrow.gllvm <- nrow(times.gllvm)

times.gllvm <- rbind(times.gllvm, tibble(gllvm=rep(NA, nrow.tot - nrow.gllvm)))
times.gllvm$gllvm[times.gllvm$gllvm > 10000] <- NA
error.gllvm <- rbind(error.gllvm, tibble(gllvm=rep(NA, nrow.tot - nrow.gllvm)))

# find number of parameters

#merge everything
times <- cbind(times, times.gllvm)
error <- cbind(error, error.gllvm)
types <- c(types, "gllvm")

rhos_small <- seq(0.01, 0.095, l=17)
rhos_large <- seq(0.1, 1, l=18)
rhos <- c(rhos_small, rhos_large)
dim <- c(rep("small", 17), rep("large", 18))




# Dplyr the tibbles for plotting

error <- error %>%
  mutate(rho = rhos, .before=1) %>%
  mutate(dim = factor(dim, levels=c("small", "large"), labels=c("small", "large")), .before=2) %>%
  pivot_longer(-c(rho, dim), names_to = "method", values_to = "error") %>%
  mutate(method=factor(method,
                       levels = types,
                       labels= types))

times <- times %>%
  mutate(rho = rhos, .before=1) %>%
  mutate(dim = factor(dim, levels=c("small", "large"), labels=c("small", "large")), .before=2) %>%
  pivot_longer(-c(rho, dim), names_to = "method", values_to = "time") %>%
  mutate(method=factor(method,
                       levels = types,
                       labels = types))

# combine the two datasets
portion_data <- cbind(times, error %>% select(error)) %>%
  pivot_longer(c("time", "error"), names_to="measure", values_to ="values")

# plot
p <- ggplot(portion_data, aes(x=rho, y=values, col=method)) +
  geom_line() +
  # scale_y_continuous(trans=log_trans(), breaks=pretty_breaks(6)) +
  scale_y_continuous(breaks=pretty_breaks(6)) +
  scale_x_continuous(breaks=pretty_breaks(20)) +
  facet_grid(rows=vars(measure), cols=vars(dim),  scale="free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid.minor = element_blank()) +
  ylab(element_blank())
  # theme(legend.position="bottom", legend.box="horizontal",
  #       strip.background = element_blank(),
  #       strip.text.y=element_blank())

p

ggsave(filename="experiments/portion_data/portion_data.png", plot=p,  width=18, height= 12, units="cm",scale=1.5)


#
#
#
#
#
#
#
#
#
#
# times <- map_dfc(types, function(type) map(bigsim, function(sim)sim[[type]]$time[1]) %>% unlist() %>% as.vector)
# colnames(times) <- types
# times <- map_dfc(times, rollmean, 3)
#
# error <- map_dfc(types, function(type) map(bigsim, function(sim)sim[[type]]$error[1]) %>% unlist() %>% as.vector)
# colnames(error) <- types
# error <- map_dfc(error, rollmean, 3)
#
#
# types_small <- c("sa", "sp", "gmf")
# types_large <- c("sa", "sp", "gmf")
#
#
#
#
#
# # small rho study
# times_small <- map_dfc(types_small, function(type) map(bigsim_small, function(sim)sim[[type]]$time[1]) %>% unlist() %>% as.vector)
# colnames(times_small) <- types_small
# # times_small <- map_dfc(times_small, rollmean, 3)
#
# error_small <- map_dfc(types_small, function(type) map(bigsim_small, function(sim)sim[[type]]$error[1]) %>% unlist() %>% as.vector)
# colnames(error_small) <- types_small
# # error_small <- map_dfc(error_small, rollmean, 3)
#
#
#
# error_small <- error_small %>%
#   mutate(rho = rho_small, .before=1) %>%
#   # REMOVING GLLVM DATA BECAUSE ERRORS ARE TOO BIG
#   mutate(gllvm=NA) %>%
#   pivot_longer(-rho, names_to = "method", values_to = "error") %>%
#   mutate(method=factor(method,
#                        levels = types_small,
#                        labels= types_small))
#
# times_small <- times_small %>%
#   mutate(rho = rho_small, .before=1) %>%
#   pivot_longer(-rho, names_to = "method", values_to = "time") %>%
#   mutate(method=factor(method,
#                        levels = types_small,
#                        labels = types_small))
#
# # combine the two datasets
# portion_data <- cbind(times_small, error_small %>% select(error)) %>%
#   pivot_longer(c("time", "error"), names_to="measure", values_to ="values")
#
# p_small <- ggplot(portion_data, aes(x=rho, y=values, col=method)) +
#   geom_line() +
#   scale_y_continuous(trans=log10_trans()) +
#   facet_grid(rows=vars(measure), scale="free_y") +
#   theme_bw()
#   # theme(legend.position="bottom", legend.box="horizontal",
#   #       strip.background = element_blank(),
#   #       strip.text.y=element_blank())
#
# p_small
#
# # large rho study
# bigsim_large <- bigsim
#
# times_large <- map_dfc(types_large, function(type) map(bigsim_large, function(sim)sim[[type]]$time[1]) %>% unlist() %>% as.vector)
# colnames(times_large) <- types_large
# # times_large <- map_dfc(times_large, rollmean, 3)
#
# error_large <- map_dfc(types_large, function(type) map(bigsim_large, function(sim)sim[[type]]$error[1]) %>% unlist() %>% as.vector)
# colnames(error_large) <- types_large
# # error_large <- map_dfc(error_large, rollmean, 3)
#
#
# error_large <- error_large %>%
#   mutate(rho = rho_large, .before=1) %>%
#   pivot_longer(-rho, names_to = "method", values_to = "error") %>%
#   mutate(method=factor(method,
#                        levels = types_large,
#                        labels= types_large))
#
# times_large <- times_large %>%
#   mutate(rho = rho_large, .before=1) %>%
#   pivot_longer(-rho, names_to = "method", values_to = "time") %>%
#   mutate(method=factor(method,
#                        levels = types_large,
#                        labels = types_large))
#
# # combine the two datasets
# portion_data <- cbind(times_large, error_large) %>%
#   pivot_longer(c("time", "error"), names_to="measure", values_to ="values")
#
# p_large <- ggplot(portion_data, aes(x=rho, y=values, col=method)) +
#   geom_line() +
#   facet_grid(rows=vars(measure), scales="free_y") +
#   theme_bw()
#
# p_large
#
# grid.arrange(p_small, p_large, nrow=1)
#
# # ALL IN 1 graph
#
# # error_small$method <- as.character(error_small$method)
# # error_large$method <- as.character(error_large$method)
# error <- rbind(error_small, error_large)
# times <- rbind(times_small, times_large)
#
# # combine the two datasets
# portion_data_all <- cbind(times, error %>% select(error)) %>% slice(1:96) %>%
#   # mutate(error_log10=log10(time)) %>%
#   # pivot_longer(c("time", "error", "error_log10"), names_to="measure", values_to ="values")
#   pivot_longer(c("time", "error"), names_to="measure", values_to ="values")
#
# p_all <- ggplot(portion_data_all, aes(x=rho, y=values, col=method)) +
#   geom_line() +
#   scale_y_continuous(trans=log10_trans())+
#   facet_wrap(vars(measure), scales="free_y", nrow = 2) +
#   xlab(TeX("$\rho$"))+
#   theme_bw() +
#   theme(axis.title.y = element_blank())
#
# p_all
#
# ggsave(filename = "./images/portion_data_error_time.png", p_all, scale=2, width=8, height=8, units="cm")

