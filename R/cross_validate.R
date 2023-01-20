#' Return lists of train and test sets to perform k
gen_kfold_sets <- function(n, k) {
  props <- (0:k)/k
  index <- props*n
  lapply(1:k, function(i){
    test <- (index[i]+1):index[i+1]
    train <- (1:n)[!(1:n %in% test)]
    list(train=train, test=test)
  })
}


crossvalidate <- function(fit, k) {
  kfold_sets <- gen_kfold_sets(fit$dimensions$n, k)
  cv_sims <- lapply(kfold_sets, function(set) {
    train <- subset(fit, set$train)
    train <- update(train, alpha=.2)

    test <- subset(fit, set$test)
    test$parameters <- train$parameter
    test <- compute_Z(test, start=test$Z, return_object = T)
    test <- compute_mean(test, return_object = T)
    test$deviance <- mean(compute_deviance(test))
    test$deviance
  })
  Reduce("+", cv_sims)/length(cv_sims)
}

if(0) {

  devtools::load_all()
  poisson  <- 0
  gaussian <- 100
  binomial <- 100
  nobs <- 1000
  q <- 4
  p <- poisson + gaussian + binomial



  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1


  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(14240)
  A <- gen_A(p, q, setting="B", prop=.4)

  fg <- gen_gllvmprime(nobs=nobs, p=p, q=q, k=k, family=family, A=A, intercept=intercept, miss.prob = 0, scale=1, phi=rep(1, p))

  fit <- gllvmprime(fg$Y, q=1, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv1 <- crossvalidate(fit, 3) # 1.58
  fit <- gllvmprime(fg$Y, q=2, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv2 <- crossvalidate(fit, 3) # 1.47
  fit <- gllvmprime(fg$Y, q=3, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv3 <- crossvalidate(fit, 3) # 1.09
  fit <- gllvmprime(fg$Y, q=4, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv4 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=5, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv5 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=6, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv6 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=7, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv7 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=8, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv8 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=9, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv9 <- crossvalidate(fit, 3) # 1.012
  fit <- gllvmprime(fg$Y, q=10, family=family, H=1, alpha=.2, maxit = 50, verbose=T)
  cv10 <- crossvalidate(fit, 3) # 1.012

  library(tidyverse)
  cv_values <- tibble(q=1:10, cv=c(1.58, 1.47, 1.09, 1.012, 1.012, 1.012, 1.012, 1.012, 1.012, 1.012))
  cv_values %>% ggplot(aes(x=q, y =cv)) +
    geom_line() +
    geom_point(size=1.6, color="white") +
    geom_point(size=1.2) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) +
    ylab(" cross-validated mean deviance") +
    theme_bw()
  ggsave("experiments/CV/CV-deviance.png")

}
