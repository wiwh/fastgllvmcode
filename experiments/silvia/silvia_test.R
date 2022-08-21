data <- read.table("experiments/silvia/dataset.dat", header=T)

F1<-factanal(x = ~x1 + x2 + x3 + x9, factors = 1, data = data[, c(1:3, 9)])

Y<-data[,c(1:3,9)]

Y<-as.matrix(Y)


Y <- as.matrix(data)
family <- "gaussian"


set.seed(3131)
fit <- fastgllvm(Y, q = 3, family="gaussian",  intercept = T, controls = list(minit=100, maxit=100,alpha=.1, eps=1e-3, learning_rate.args=list(method="constant")))
fit <- update(fit)

fit2 <- ffa(Y, 3, iteratively_update_Psi = T)

plot(fit2$A, psych::Procrustes(fit$parameters$A, fit2$A)$loadings); abline(0,1,col=2)
plot(fit2$Psi, fit$parameters$phi); abline(0,1,col=2)
fit$parameters$A
F1$loadings
