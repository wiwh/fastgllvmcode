library(tidyverse)
ARE <- readRDS("experiments/ARE/ARE_n_1e4_seed1231.rds")
ARE <- sapply(1:4, function(i) unlist(ARE[i,]))
colnames(ARE) <- c("PROMES", "PROMES2", "S-PROMES", "SPROMES2")
ARE <- as_tibble(ARE) %>% mutate(p=c(5, 10, 20, 30, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500)) %>%
  pivot_longer(-p, names_to="Estimator", values_to = "ARE")


ggplot(ARE %>% filter(Estimator %in% c("PROMES", "S-PROMES")), aes(x=p,y=ARE, col=Estimator, shape=Estimator))+
  geom_line() +
  geom_point(size=2, col="white") +
  geom_point(size=1.5) +
  scale_y_continuous(limits = c(.5,1)) +
  scale_x_continuous(breaks=unique(ARE$p), labels = unique(ARE$p), limits=c(5, 500), minor_breaks = NULL, trans="sqrt") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45))


ggsave("experiments/ARE/ARE_plot.png", width = 12, height = 7)
