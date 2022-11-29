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
    train <- update(train)

    test <- subset(fit, set$test)
    test$parameters <- train$parameters
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
q <- 2
p <- poisson + gaussian + binomial



intercept <- T
k <- 1
if(k==0 & intercept) k <- 1


family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
set.seed(14240)
fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1, phi=rep(1, p))

fit <- fastgllvm(fg$Y, q=1, family=family, H=1, alpha=.2, maxit = 50)
cv1 <- crossvalidate(fit, 3)

fit <- fastgllvm(fg$Y, q=2, family=family, H=1, alpha=.2, maxit = 50)
cv2 <- crossvalidate(fit, 3)
fit <- fastgllvm(fg$Y, q=3, family=family, H=1, alpha=.2, maxit = 50)
cv3 <- crossvalidate(fit, 3)
fit <- fastgllvm(fg$Y, q=4, family=family, H=1, alpha=.2, maxit = 50)
cv4 <- crossvalidate(fit, 3)
fit <- fastgllvm(fg$Y, q=5, family=family, H=1, alpha=.2, maxit = 50)
cv5 <- crossvalidate(fit, 3)


plot(fit)


}
