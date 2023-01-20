library(tidyverse)
# this are obtained as gglvm where first 10 are poisson, next 10 are binary, loadings in seq(-2, 2), with a bias in loading of +.5

# Now gradient graphs
sims_full <- readRDS(file="./experiments/gradients/sims_gradient_full.rds")
sims_simple <- readRDS(file="./experiments/gradients/sims_gradient_simple.rds")


prepare_data <- function(gradient_sims) {
  gradient_sims_poisson <- as_tibble(gradient_sims) %>% select(1:10) %>% pivot_longer(everything(), names_to = "loading", values_to="value") %>% mutate(family="Poisson")
  gradient_sims_binary  <- as_tibble(gradient_sims) %>% select(11:20) %>% pivot_longer(everything(), names_to = "loading", values_to="value") %>% mutate(family="Bernoulli")
  gradient_sims <- rbind(gradient_sims_poisson, gradient_sims_binary) %>%
    mutate(family=factor(family, labels=c("Poisson", "Bernoulli"), levels=c("Poisson", "Bernoulli")))%>%
    mutate(loading = factor(loading, levels = paste0("V", 1:20), labels=paste0("", 1:20)))

  gradient_sims
}

plot_sim <- function(gradient_sims, gradient_sims_rescaled, scale_y_limits=c(-2, 10)) {
  gradient_sims <- prepare_data(gradient_sims) %>% mutate(rescaled="Original")
  gradient_sims_rescaled <- prepare_data(gradient_sims_rescaled) %>% mutate(rescaled="Re-scaled")

  gradient_sims <- rbind(gradient_sims, gradient_sims_rescaled)

  gradient_sims %>% ggplot(aes(x=loading, y=value)) +
    geom_boxplot(outlier.shape=NA, size=.2) +
    scale_y_continuous(limits=scale_y_limits, breaks=c(0,1,4,8,12)) +
    facet_grid(rescaled ~family, scale="free_x") +
    geom_hline(yintercept=0, colour="red", size=.2) +
    geom_hline(yintercept=1, colour="red", linetype="dashed", size=.2) +
    theme_bw() +
    ylab("") +
    xlab("Loading")
}

plot_sim(-sims_full$sim_full$AB, sims_full$sim_full_rescaled$AB,  scale_y_limits = c(-2, 4))
ggsave(filename="./experiments/gradients/sims_gradient_full.png", width = 14*.7, height=8*.7)
plot_sim(-sims_simple$sim_simple$AB, sims_simple$sim_simple_rescaled$AB,  scale_y_limits = c(-2, 12))
ggsave(filename="./experiments/gradients/sims_gradient_simple.png", width=14*.7, height =8*.7)
