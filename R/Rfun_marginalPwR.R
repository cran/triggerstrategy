# Rfun_marginalPwR
# 2017-12-08
#
#' @name marginalPwR
#' @title Marginal Power Rate
#' @description This function computes the marginal powers.
#' @param cvec a vector of critical boundaries
#' @param t a vector of information times
#' @param delta a number shows the drift paramter
#' @return a number shows the marginal power (delta isn't equal to zero) or type I error (delta is zero)
#' @export
#' @import mvtnorm
#' @author Jiangtao Gou
#' @examples
#' marginalPwR(c(2.218,2.218),c(0.1,0.5,1),delta=3)
#' marginalPwR(1.96,t=1,delta=3)
#
marginalPwR <- function (cvec,t,delta=0) {
  K <- min(length(cvec),length(t))
  if ( K <= 1) {
    pwr <- pnorm(q=cvec[1], mean=delta*sqrt(t[1]), sd=1, lower.tail=FALSE)
    return(pwr)
  } # Add this if condition on 2020-03-05 08:10
  # tvec <- t[1:K]/t[K] # Don't standardize t!!!
  tvec <- t[1:K]
  lowerB <- rep(-Inf,times=K)
  upperB <- cvec[1:K]
  meanV <- delta*sqrt(tvec)
  corrM <- corrMatGenerator(tp=tvec,ts=vector(mode="numeric",length=0),rhops=1)
  resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
  result <- resultIntgl[1]
  pwr <- 1-result
  return (pwr)
}#
