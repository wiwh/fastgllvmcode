library(tidyverse)
source("./experiments/simu_functions.R")
load("./experiments/low_dimensions/low_dimensions_results.Rdata")
# re-organize the results from setting > iterations > methods to setting > method > iterations


# methods:
methods  <- c("sa", "sp","gllvm", "gmf")

# Example: for simulation number i
i <- 1
gllvm <- extract_loadings(bigsim, i, "gllvm")
gllvm.time <- extract_time(bigsim, i, "gllvm")

sa <- extract_loadings(bigsim, i, "sa")
sa.time <- extract_time(bigsim, i, "sa")

sp <- extract_loadings(bigsim, i, "sp")
sp.time <- extract_time(bigsim, i, "sp")

gmf <- extract_loadings(bigsim, i, "gmf")
gmf.time <- extract_time(bigsim, i, "gmf")

p <- bigsim[[i]][[i]]$p
true_loadings <- c(seq(2,-2, l=p), seq(-1, 1, l=p))
plot(true_loadings, colMeans(sa)-true_loadings, )
points(true_loadings, colMeans(sp)-true_loadings)
points(true_loadings, colMeans(gmf)-true_loadings)


# draw boxplots

boxplot(sa, outline=T, main=sa.time)
points(true_loadings, col=2)

boxplot(sp, outline=T, main=sp.time)
points(true_loadings, col=2)

boxplot(gmf, outline=T, main=gmf.time)
points(true_loadings, col=2)

boxplot(gllvm, outline=T, main=gllvm.time)
points(true_loadings, col=2)

# Produce graph outputs


# Produce table outputs

# Get the table of median of the mean MSE
mse.median <- t(sapply(1:6, function(i){
  sapply(methods, function(est){
    loadings <- extract_loadings(bigsim, i, est)
    mse <- median(apply(loadings, 1, function(est)sqrt(mean((true_loadings - est )^2))))
    round(mse, 2)
  })
}))

mse.mean <- t(sapply(1:6, function(i){
  sapply(methods, function(est){
    loadings <- extract_loadings(bigsim, i, est)
    mse <- mean(apply(loadings, 1, function(est)sqrt(mean((true_loadings - est )^2))))
    round(mse,2)
  })
}))

# combine the two

table.mse <- sapply(seq_along(methods), function(i){
  apply(cbind(mse.mean[,i], mse.median[,i]),1, function(res)paste0(res[1]," (",res[2], ")" ))
})

settings <- extract_settings(bigsim)

table.mse <- cbind("", settings[,2], table.mse)

(table.mse.tex <- knitr::kable(table.mse, format="latex", booktabs = T))


# get the table of all timees
table.times <- t(sapply(1:6, function(i){
  round(sapply(methods, function(est) extract_time(bigsim, i, est)),2)
}))

table.times <- cbind("", settings[,2], table.times)
(table.times.tex <- knitr::kable(table.times, format="latex", booktabs=T))




