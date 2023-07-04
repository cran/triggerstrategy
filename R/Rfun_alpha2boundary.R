# Rfun_alpha2boundary
# 2020-03-02 13:20-14:20
#' @name alpha2boundary
#' @title Convert cumulative alpha levels to normal critical boundaries
#' @description This function converts cumulative alpha levels into normal critical boundaries.
#' @param alphas a list of cumulative errors from some error spending functions
#' @param t a vector of information times
#' @param initIntvl a pair of numbers as the lower and upper bounds of critical boundaries, used for \code{stats::uniroot} function
#' @return a vector of critical boundaries
#' @export
#' @import stats
#' @import mvtnorm
#' @import ldbounds
#' @author Jiangtao Gou
#' @details The current version of \code{ldbounds::ldBounds} may not work for Hwang-Shih-DeCani boundaries.
#' @examples
#' library(ldbounds)
#' tvec <- c(0.5,1)
#' result <- ldbounds::ldBounds(t=tvec, iuse=1, alpha=0.05, sides=1)
#' print(result$upper.bounds)
#' bd <- alpha2boundary(alphas = result$exit.pr, t=tvec)
#' print(bd)
#
alpha2boundary <- function(alphas, t, initIntvl=c(1,2*stats::qnorm(p=alphas[1], lower.tail = FALSE))) {
  K <- length(alphas)
  Kt <- length(t)
  bounds <- rep(0, times=K)
  if (K != Kt | K < 2) {
    return(bounds)
  }
  bounds[1] <- stats::qnorm(p=alphas[1],lower.tail = FALSE)
  #
  target <- function(newbound, currentbound, alphasub, tsub) {
    Ksub <- length(tsub)
    lowerB <- rep(-Inf,times=Ksub)
    upperB <- c(currentbound, newbound)
    meanV <- rep(0,Ksub)
    corrM <- corrMatGenerator(tp=tsub,ts=vector(mode="numeric",length=0),rhops=1)
    resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
    result <- resultIntgl[1]
    return(result - 1 + alphasub[Ksub])
  }
  #
  for (stage in 2:K) {
    result <- stats::uniroot(target, lower=initIntvl[1], upper=initIntvl[2], tol=2.5e-16, currentbound=bounds[1:(stage-1)], alphasub = alphas[1:stage], tsub = t[1:stage])
    bounds[stage] <- result$root
  }
  return(bounds)
}#
