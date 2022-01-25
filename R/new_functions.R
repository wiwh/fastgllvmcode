#' Function factory to return a function that generates a learning rate as a function of i
#'
#' @param method: the name of the method to use. One of "constant", "lin", "exp", "linexp", "spall".
#' @param method.args: list of method arguments

ff_learning_rate <- function(method="constant", maxit, learning_rate.args){
  if(is.null(learning_rate.args)) learning_rate.args <- list()
  if(is.null(learning_rate.args$start)) learning_rate.args$start <- 40
  if(is.null(learning_rate.args$end)) learning_rate.args$end <- 0.1
  if(is.null(learning_rate.args$constant)) learning_rate.args$constant <- 10

  if(method=="constant"){
    learning_rate <- function(i){
      rep(learning_rate.args$constant, length(i))
    }
  }
  if(method=="lin"){
    lr <- with(learning_rate.args,
               seq(start, end, l=maxit))
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(method=="exp"){
    lr <- with(learning_rate.args, {
      rho <- exp(log(end/start)/(maxit-1))
      rep(start, maxit) * rho**(0:(maxit-1))
    })
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(method=="linexp"){
    if(is.null(learning_rate.args$weight.exp))learning_rate.args$weight.exp <- .5
    lr.lin <- with(learning_rate.args,
               seq(start, end, l=maxit))
    lr.exp <- with(learning_rate.args, {
      rho <- exp(log(end/start)/(maxit-1))
      rep(start, maxit) * rho**(0:(maxit-1))
    })
    lr <- (1-learning_rate.args$weight.exp) * lr.lin + learning_rate.args$weight.exp * lr.exp
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(method=="spall"){
    if(is.null(learning_rate.args$rate)) learning_rate.args$rate <- 10
    lr <- with(learning_rate.args, {
      b <- exp(log(end / start) / rate)
      b <- (maxit * b - 1) / (1 - b)
      start*((1+b)/(1:maxit + b))**rate
    })
    learning_rate <- function(i){
      lr[i]
    }
  }
  learning_rate
}

#
# maxit <- 1000
# learning_rate <- ff_learning_rate(method="spall", learning_rate.args = list(start=40, end=0.1,rate=10))
# plot(learning_rate(1:maxit), type="l")
#
# learning_rate <- ff_learning_rate(method="linexp", learning_rate.args = list(start=40, end=0.1, weight.exp=.9))
# lines(learning_rate(1:maxit), col=2)
#
# learning_rate <- ff_learning_rate(method="exp", learning_rate.args = list(start=40, end=0.1))
# lines(learning_rate(1:maxit), col=3)
