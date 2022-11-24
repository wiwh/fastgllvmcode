library(tidyverse)
library(scales)
library(gridExtra)
library(zoo)
library(latex2exp)

if(0){
  # here we want to construct a data-set with columns loading_id, loading_true, loading_est, n, sim_num,
  # with all observations as rows

  simnames <- paste0("./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_",
                     c("1-20",
                       "21-40",
                       "41-60",
                       "61-80",
                       "81-100"),
                     ".Rdata")

  # for each of the sim repetitions (20 per file), there are 10 settings of n, for each of them, we want to extract the loading and the true
  # value of the loadings

  sim <- lapply(seq_along(simnames), function(i){
    load(simnames[i])
    sim <- lapply(seq_along(bigsim), function(j){
      sim <- map(bigsim[[j]], function(b){
        load.true <- as.vector(b$dat.par$A)
        load.est  <- as.vector(b$sa$A)
        tibble(sim.id=(i-1)*20+j,  sim.n=b$dat.n, load.id=1:length(load.true), load.true=load.true, load.est=load.est, time=b$sa$time[1])
      })
      do.call(rbind, sim)
    })
    do.call(rbind, sim)
  })
  sim <- do.call(rbind, sim)
  save(sim, file="experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_tidy_new.Rdata")
}

load("experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_tidy_new.Rdata")

# compute bias and variance
mse_bias_var <- sim %>%
  group_by(sim.n, load.true) %>%
  summarise(mse=(mean((load.true - load.est)^2)), bias=abs(mean(load.est) - mean(load.true)), sd=sd(load.est))

bias <- ggplot(mse_bias_var, aes(sim.n)) +
  geom_line(aes(y=bias)) +
  geom_line(aes(y=sd)) +
  facet_grid(load.true ~ ., scale="free_y") +
  scale_x_continuous(breaks=pretty_breaks(10))

mse_bias_var_avg <- mse_bias_var %>% group_by(sim.n) %>% summarise_all(mean)
mse_bias_var_avg <- mse_bias_var %>%
  subset(load.true=="1") %>%
  select(-mse) %>%
  pivot_longer(c(bias, sd), names_to="measure", values_to = "value")

bias_avg <- ggplot(mse_bias_var_avg, aes(sim.n, value, col=measure)) +
  geom_line() +
  xlab("n") +
  ylab("")+
  geom_point()+
  scale_x_continuous(breaks=pretty_breaks(10)) +
  scale_color_discrete(labels=c("Bias", "Variance"))+
  theme_bw() +
  theme(legend.position = c(.9, 0.8),
        legend.title = element_blank())

bias_avg

sim$load.true <- factor(sim$load.true, levels = sort(unique(sim$load.true), decreasing=T), labels=as.character(sort(unique(sim$load.true), decreasing=T)))

#now the boxplots

boxplots <- ggplot(sim, aes(sim.n,  load.est)) +
  geom_blank() +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(load.true ~., scales="free_y") +
  geom_boxplot(aes(group=sim.n, colour=load.true)) +
  scale_x_continuous(breaks=pretty_breaks(10))+
  scale_y_continuous(breaks=pretty_breaks(1)) +
  xlab("n")+
  ylab("loading value") +
  theme_bw() +
  theme(panel.spacing = unit(0, "mm"),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

boxplots

bias_avg


grid.arrange(boxplots, bias_avg, nrow=2)

ggsave(filename="experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_boxplots.png", plot=boxplots,  width=18, height= 16, units="cm",scale=1.5)
ggsave(filename="experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_bias_var.png", plot=bias_avg,  width=18, height= 12, units="cm",scale=1.5)
