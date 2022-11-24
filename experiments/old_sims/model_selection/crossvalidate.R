library(caret)
#https://remiller1450.github.io/s230f19/caret2.html

n <- 1000
p <- 100
q <- 10 # number of true variables
# gen data
X <- matrix(rnorm(p*n), n, p)
colnames(X) <- paste0("X", 1:p)
b <- c(rnorm(q), rep(0, p-q))
linpar <- X %*% b
family <- binomial(link="logit")

Y <- rbinom(n, 1, family$linkinv(linpar))

# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(x=X, y=Y, trControl=train_control, method="glm", family="binomial")
# summarize results
print(model)
