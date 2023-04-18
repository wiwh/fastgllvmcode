library(tidyverse)

# Generate a random sparse loading matrix
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


plot_sparse_structure <- function(A) {
  colnames(A) <- 1:ncol(A)

  A.tbl <- as_tibble(A) %>%
    mutate(y=1:nrow(A)) %>%  # create the y axis
    pivot_longer(-y, names_to="x", values_to="z") %>%
    mutate(x=as.numeric(x)) %>%
    mutate(z=ifelse(z!=0, 1, 0)) %>%  # nonzero loadings are transformed to 1
    mutate(z=factor(z, labels=c("loading is non-zero","loading is zero"), levels=1:0))

  p <- A.tbl %>% ggplot(aes(x=x, y=y, fill=z)) +
    geom_raster() +
    scale_y_reverse() +
    scale_fill_brewer(palette = "Greys", direction=-1)+
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(legend.position="none")+
    xlab("Latent Variable") +
    ylab("Response") +
    theme(text=element_text(size=18))
  p
}


# Example
# -------

# Generate the loading matrix (can be replaced by any)
A <- gen_A(500, 5, nonzero=200)

# Plot it
p <- plot_sparse_structure(A)
p
