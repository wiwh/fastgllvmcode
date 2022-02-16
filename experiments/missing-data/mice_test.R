devtools::load_all()
library(missMDA)
library(FactoMineR)

data(ozone)
X <- ozone[,1:11]
M <- is.na(dat) *1
X <- imputePCA(X)
X.pca <- PCA(X, ncp=3)$ind$coord



fit <- fastgllvm(M, X=cbind(1, X.pca), q=3, learning_rate.args = list(start=2, end=0.5))

Ysim <- simulate_fastgllvm(fit, n=10000)
plot(cov(M), cov(Ysim$Y))




