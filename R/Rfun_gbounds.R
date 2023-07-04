# Rfun_gbounds
# 2020-03-04
# 2020-03-21 Add original POC and OBF
# 2021-06-21 Programming HSD directly instead of using R package gsDesign
#
#' @name gbounds
#' @title Critical boundary in group sequential trials
#' @description This function computes the critical boundaries and the error spent until each stage in group sequential trials
#' @param t a vector of information times
#' @param iuse a number of the type of the error spending function, from -2, -1, 1, 2, 3, 4
#' @param alpha a number of type I error rate
#' @param phi a parameter for the power family or the HSD gamma family
#' @return a list of two vectors: \code{bd} critical boundaries, \code{er} error spent until each stage
#' @export
#' @import stats
#' @import ldbounds
#' @import mvtnorm
#' @author Jiangtao Gou
#' @details If the original Pocock is implemented, we specify \code{iuse=-2}. If the original OBrien-Flemming is implemented, we specify \code{iuse=-1}.
#' @examples
#' t<-c(0.5,0.8,1)
#' iuse <- 4
#' gbounds(t=t, iuse=iuse)
#' gbounds(t=(1:5)/5, iuse=4, alpha=0.01, phi=-4)
#' gbounds(t=(1:5)/5, iuse=-2, alpha=0.01)
#
gbounds <- function(t, iuse=1, alpha=0.05, phi=rep(1,length(alpha))) {
  t <- t/t[length(t)] # Normalization
  if (length(t) <= 1) {
    er <- alpha
    bd <- stats::qnorm(p=alpha,lower.tail = FALSE)
    result <- list(bd=bd, er=er)
    return(result)
  }
  if (iuse == 1 | iuse == 2 | iuse == 3) {
    ldresult <- ldbounds::ldBounds(t=t, iuse=iuse, alpha=alpha, phi=phi, sides=1)
    bd <- ldresult$upper.bounds
    er <- ldresult$exit.pr
    result <- list(bd=bd, er=er)
    return(result)
  } else if (iuse == 4) {
    # sfresult <- gsDesign::sfHSD(alpha=alpha, t=t, param=phi)
    # er <- sfresult$spend
    #
    er <- if (phi == 0) t * alpha else alpha * (1. - exp(-t * phi)) / (1 - exp(-phi))
    #
    bd <- alpha2boundary(alphas=er, t=t)
    result <- list(bd=bd, er=er)
    return(result)
  } else if (iuse == -2) {
    # Orginal Pocock
    target <- function (c,t,alpha) {
      K <- length(t)
      lowerB <- rep(-Inf,K)
      upperB <- c
      meanV <- rep(0,K)
      corrM <- corrMatGenerator(tp=t,ts=vector(mode="numeric",length=0),rhops=1)
      result <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
      return(result - 1 + alpha)
    }
    rootresult <- stats::uniroot(target, lower = qnorm(p=alpha, lower.tail=FALSE), upper = 10, tol = 2.5e-16, t=t, alpha = alpha)
    bd <- rep(rootresult$root, times=length(t))
    er <- boundary2alpha(cvec=bd, t=t)
    result <- list(bd=bd, er=er)
    return(result)
  } else if (iuse == -1) {
    # Orginal OBF (1979)
    target <- function (c,t,alpha) {
      K <- length(t)
      lowerB <- rep(-Inf,K)
      upperB <- c/sqrt(t)
      meanV <- rep(0,K)
      corrM <- corrMatGenerator(tp=t,ts=vector(mode="numeric",length=0),rhops=1)
      result <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
      return(result - 1 + alpha)
    }
    rootresult <- stats::uniroot(target, lower = qnorm(p=alpha, lower.tail=FALSE), upper = 10, tol = 2.5e-16, t=t, alpha = alpha)
    bd <- rootresult$root/sqrt(t)
    er <- boundary2alpha(cvec=bd, t=t)
    result <- list(bd=bd, er=er)
    return(result)
  } else {
    result <- list(bd=rep(0,times=length(t)), er=rep(0,times=length(t)))
    return(result)
  }
}#
