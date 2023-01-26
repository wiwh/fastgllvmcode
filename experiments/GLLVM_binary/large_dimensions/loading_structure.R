# gen A according to the two sparsity settings
library(tidyverse)
library(gridExtra)
library(ggpubr)


gen_A <- function(p, q, setting="A", nonzero=100, prop=.4) {
  # setting "A" has only the 100 first loadings that are non-zero
  A <- matrix(runif(p*q, -2, 2), p, q)
  if (setting == "A") {
    if(p>nonzero) {
      A[(nonzero+1):p, ] <- 0
    }
  } else if (setting =="B") {
    shift_size = round(p * (1-prop)/(q-1))

    for (k in 1:q) {

      nonzero_start = (k-1) * shift_size + 1

      if (k == q) {
        nonzero_end = p
      } else {
        nonzero_end   = nonzero_start + round(prop*p) - 1
      }

      nonzeros <- (1:p) %in% (nonzero_start:nonzero_end)
      A[!nonzeros,k] <- 0
    }
  }
  A
}




gen_graph_data <- function(q, p, prop) {
  AB.tbl <- lapply(p, function(p) {
    # Setting A
    A <- gen_A(p, q, "A", prop=prop)
    colnames(A) <- 1:q

    A.tbl <- as_tibble(A) %>% mutate(y=1:nrow(A)) %>% pivot_longer(-y, names_to="x", values_to="z") %>% mutate(x=as.numeric(x)) %>% mutate(z=ifelse(z!=0, 1, 0)) %>% mutate(Setting="Setting A")
    #mutate(x=factor(x, levels=1:5, labels=1:5), y=factor(y, levels=1:p, labels=1:p))

    # Setting B
    A <- gen_A(p, q, "B", prop=prop)
    colnames(A) <- 1:q
    B.tbl <- as_tibble(A) %>% mutate(y=1:nrow(A)) %>% pivot_longer(-y, names_to="x", values_to="z") %>% mutate(x=as.numeric(x)) %>% mutate(z=ifelse(z!=0, 1, 0)) %>% mutate(Setting="Setting B")

    AB.tbl <- rbind(A.tbl, B.tbl) %>% mutate(p= paste0("p = ", p))

  })
  AB.tbl <- do.call(rbind, AB.tbl)
}

Ab.tbl <- gen_graph_data(5, c(100, 200, 300), .4)

pA <- AB.tbl %>% filter(Setting=="Setting A") %>% mutate(z=factor(z, labels=c("loading is non-zero","loading is zero"), levels=1:0)) %>%
  ggplot(aes(x=x, y=y, fill=z)) +
  geom_raster() +
  scale_y_reverse() +
  facet_grid(.~ p) +
  scale_fill_brewer(palette = "Greys", direction=-1)+
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.position="none")+
  xlab("k") +
  ylab("j") +
  theme(text=element_text(size=18)) +
  ggtitle("Setting A")

pB <- AB.tbl %>% filter(Setting=="Setting B") %>% mutate(z=factor(z, labels=c("loading is non-zero","loading is zero"), levels=1:0)) %>%
  ggplot(aes(x=x, y=y, fill=z)) +
  geom_raster() +
  scale_y_reverse() +
  facet_grid(.~ p) +
  scale_fill_brewer(palette = "Greys", direction=-1)+
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(.2, .2)) +
  xlab("k") +
  ylab("j") +
  theme(text=element_text(size=18)) +
  ggtitle("Setting B")

p <- ggarrange(pA, pB, nrow=1, common.legend = TRUE, legend="bottom")
ggsave("experiments/GLLVM_binary/large_dimensions/loading_structure.png", plot=p,  height=8, width=12)

# this is for chapter 2
Ab.tbl.2 <- gen_graph_data(4, 200, .5)
pC <- Ab.tbl.2 %>% filter(Setting=="Setting B") %>% mutate(z=factor(z, labels=c("loading is non-zero","loading is zero"), levels=1:0)) %>%
  ggplot(aes(x=x, y=y, fill=z)) +
  geom_raster() +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))+
  scale_fill_brewer(palette = "Greys", direction=-1)+
  theme_bw()  +
  xlab("q") +
  ylab("loading") +
  theme(text=element_text(size=40),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.title=element_blank())
pC

ggsave("experiments/CV/loading_structure_chap2.png", height=14, width=10)


