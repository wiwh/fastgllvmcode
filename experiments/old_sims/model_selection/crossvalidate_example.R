library(tidyverse)
library(ggplot2) ## For plotting
library(caret)   ## For model fitting and evaluation
library(visreg)  ## For visualizing regression models
library(plotROC) ## For constructing ROC curves
library(mgcv)    ## For fitting GAM models
library(kernlab) ## Contains an example dataset
library(glmnet)  ## For fitting regularized models

## Load the Donner Party data and plot it

n <- 1000
p <- 100
q <- 10 # number of true variables
# gen data
X <- matrix(rnorm(p*n), n, p)
colnames(X) <- paste0("X", 1:p)
b <- c(rnorm(q), rep(0, p-q))
linpar <- X %*% b
family <- binomial(link="logit")

Y <- rbinom(n, 1, family$linkinv(linpar)) %>% as.factor()
dat <- as_tibble(cbind(Y, X))

dat$Y <- as.factor(dat$Y)
fit.control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

set.seed(123)
fit <- train(y=Y, x=X[,1:(q)], method = "glm",
             family = "binomial", trControl = fit.control)

fit
