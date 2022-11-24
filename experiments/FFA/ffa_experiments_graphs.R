library(tidyverse)
readRDS("./experiments/FFA/sims_speed_fad_3.rds")

# sim setting > sim rep > sim results

extract_info <- function(sims) {
  simres <- lapply(sims, function(sim) t(sapply(sim, function(sim_i) {
    c(n=sim_i$n, p=sim_i$p, q=sim_i$q, r=sim_i$r, time.fad=sim_i$fad.time[3], time.ffa=sim_i$ffa.time[3], error.fad=sim_i$fad.error.equality, error.ffa=sim_i$ffa.error.equality)
  })))

  tbl <- as_tibble(do.call(rbind, simres))
  colnames(tbl) <- c("n", "p", "q", "r", "time.fad.elapsed","time.ffa.elapsed","error.fad","error.ffa")
  tbl
}


r3 <- readRDS("./experiments/FFA/sims_speed_fad_3.rds") %>% extract_info()
r5 <- readRDS("./experiments/FFA/sims_speed_fad_5.rds") %>% extract_info()

s <- rbind(r3, r5)

s <- s %>% mutate(npq= paste0("n = ", n, ", p = ", p)) %>%
  mutate(speedup=time.fad.elapsed/time.ffa.elapsed) %>%
  select(-n, -p, -time.fad.elapsed, -time.ffa.elapsed) %>%
  mutate(q=paste0("q = ", q))

ggplot(s, aes(x=r, y=speedup, col=npq, group=r)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(q~npq, scales = "free") +
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(size=14)) +
  ylab("Speed-up")+
  scale_y_continuous(breaks=c(1, 10, 20, 40, 40, 50, 100, 150), minor_breaks = NULL) +
  geom_hline(yintercept=1) +
  scale_x_continuous(breaks=1:10)


ggsave("./experiments/FFA/ffa_speedup.png")

