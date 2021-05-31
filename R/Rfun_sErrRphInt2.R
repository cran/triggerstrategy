# Rfun_sErrRphInt2
# Partial hierarchical, error rate of the intersection hypothesis, two endpoints
# NOTE: sErrRphInt written on 2017-12-08 has typos and can't work correctly when Ks - Kp >=2. 
# 2020-03-04
#
#' @name sErrRphInt2
#' @title Type I error rate of the overall null hypothesis using the partial hierarchical design
#' @description This function computes the type I error rate of the overall null hypothesis using the partial hierarchical group sequential design.
#' @param cvec0 a vector of critical boundaries for testing H0
#' @param cvec1 a vector of critical boundaries for testing H1
#' @param t0 a vector of information times for H0
#' @param t1 a vector of information times for H1
#' @param tc0 a vector of calendar times for H0
#' @param tc1 a vector of calendar times for H1
#' @param rho a value of the correlation between the test statistics for H0 and H1. 
#' @return a number shows the type I error rate of testing H0 intersect H1
#' @export
#' @import mvtnorm
#' @author Jiangtao Gou
#' @examples
#' alpha0 <- 0.03
#' alpha1 <- 0.02
#' iuse0 <- 4
#' iuse1 <- 4
#' phi0 <- -4
#' phi1 <- 1
#' tc0 <- c(3,6,9,12)
#' tc1 <- c(6,12,18,24)
#' t0 <- c(0.3,0.6,0.9,1)
#' t1 <- (1:4)/4
#' rho <- 0
#' cvecList0 <- gbounds(t=t0,iuse=iuse0, 
#'     alpha=alpha0,phi=phi0)
#' cvec0 <- cvecList0$bd
#' cvecList1 <- gbounds(t=t1,iuse=iuse1, 
#'     alpha=alpha1,phi=phi1)
#' cvec1 <- cvecList1$bd
#' result <- sErrRphInt2(cvec0, cvec1, 
#'     t0, t1, tc0, tc1, rho)
#' print(result)
#
sErrRphInt2 <- function (cvec0, cvec1, t0, t1, tc0=t0, tc1=t1, rho=0) {
  stageK0 <- length(t0)
  stageK1 <- length(t1)
  # Marginal type I error for rejecting H0
  alpha0 <- marginalPwR(cvec=cvec0,t=t0,delta=0)
  # Testing H0 lasts longer than testing H1
  if (is.na(which(tc1 >= tc0[stageK0])[1])) {
    #
    meanV <- rep(0,times=stageK0+1)
    corrM <- corrMatGenerator(tp=t0, ts=t1[stageK1],rhops=rho)
    lowerB <- c(rep(-Inf,times=stageK0),cvec1[stageK1])
    upperB <- c(cvec0,Inf)
    tempIntgl <- mvtnorm::pmvnorm(lowerB,upperB,meanV,corrM,algorithm=Miwa(steps=128))
    alphaDiff <- tempIntgl[1]
    alpha01 <- alpha0 + alphaDiff
    return(alpha01)
  }
  # Testing H1 lasts longer than testing H0
  idx_t1_tail <- which(tc1 >= tc0[stageK0]) # Pick the stages of H1 which are not skipped
  stageK1_tail <- length(idx_t1_tail)
  t1_tail <- t1[idx_t1_tail]
  cvec1_tail <- cvec1[idx_t1_tail]
  #
  tempResult <- rep(0,times=stageK1_tail)
  #
  for (i in 1:stageK1_tail) {
    meanV <- rep(0,times=stageK0 + i)
    corrM <- corrMatGenerator(tp=t0, ts=t1_tail[1:i], rhops=rho)
    if (i == 1) {
      lowerB <- c(rep(-Inf,times=stageK0), cvec1_tail[i])
      upperB <- c(cvec0,Inf)
    } else {
      lowerB <- c(rep(-Inf,times=stageK0+i-1), cvec1_tail[i])
      upperB <- c(cvec0, cvec1_tail[1:(i-1)], Inf) 
    }
    tempIntgl <- mvtnorm::pmvnorm(lowerB,upperB,meanV,corrM,algorithm=Miwa(steps=128))
    tempResult[i] <- tempIntgl[1]
  }
  alphaDiff <- sum(tempResult)
  alpha01 <- alpha0 + alphaDiff
  return(alpha01)
}# 
