library(tidyverse)
dt  <- readRDS("experiments/FFA/likelihood_improvement.rds")

dt <- dt %>% mutate(i=2:(n()+1))
dt <- rbind(c(-730438, -730438, 1), dt)

dt <- dt %>% pivot_longer(-i, values_to = "likelihood", names_to="Algorithm")


ggplot(dt, aes(x=i, y=likelihood, col=Algorithm)) +
  geom_line() +
  xlab("Iteration") +
  ylab("") +
  geom_hline(yintercept=max(dt$likelihood), lty="dotted") +
  scale_x_continuous(breaks=c(1,50,100,150,200,250,300,350,400,450,500), minor_breaks = NULL) +
  scale_y_continuous(breaks=c(max(dt$likelihood)), label="Maximum \n log-likelihood") +
  # geom_text(x = 20, y=max(dt$likelihood)+200, label=paste0("Max. loglikelihood = ", max(dt$likelihood))) +
  theme_bw(base_size=20)

ggsave("experiments/FFA/ffa_likelihood_improvements.png", width=11, height=5)

