# Rfun_sPwRph2.R
# 2020-03-05
# A deep extension of sPwRph.R, version 2017-12-08
#
#' @name sPwRph2
#' @title Power of testing the secondary hypothesis in partially hierarchical design
#' @description This function computes the power of testing the secondary hypothesis in partially hierarchical design.
#' @param cvec0 a vector of critical boundaries for testing H0
#' @param cvec1 a vector of critical boundaries for testing H1
#' @param delta0 a value of drift parameter for testing H0
#' @param delta1 a value of drift parameter for testing H1
#' @param t0 a vector of information times for H0
#' @param t1 a vector of information times for H1
#' @param tc0 a vector of calendar times for H0
#' @param tc1 a vector of calendar times for H1
#' @param rho a value of correlation coefficient between H0 and H1
#' @return a value of the probability that H1 is rejected, the power
#' @export
#' @examples 
#' alpha <- 0.05
#' alpha0 <- 0.03
#' iuse0 <- 4
#' iuse1h <- 4
#' iuse1t <- 4
#' phi0 <- -4
#' phi1h <- 1
#' phi1t <- 1
#' tc0 <- c(3,6,9,12)
#' tc1 <- c(6,12,18,24)
#' t0 <- c(0.3,0.6,0.9,1)
#' t1 <- (1:4)/4
#' rho <- 0
#' cvecList0 <- gbounds(t=t0, iuse=iuse0, 
#'     alpha=alpha0, phi=phi0)
#' cvec0 <- cvecList0$bd
#' cvecList1 <- sBoundsPh2(alpha, alpha0, 
#'     t0, t1, tc0, tc1, 
#'     rho, iuse0, iuse1h, iuse1t, 
#'     phi0, phi1h, phi1t)
#' cvec1 <- cvecList1$bd
#' sPwRph2(cvec0, cvec1, 
#'     delta0=2, delta1=3, 
#'     t0, t1, tc0, tc1, 
#'     rho=0) 
#
sPwRph2 <- function(cvec0, cvec1, delta0, delta1, t0, t1, tc0=t0, tc1=t1, rho=0) {
  stageK0 <- length(t0)
  stageK1 <- length(t1)
  pw1max <- marginalPwR(cvec=cvec1, t=t1, delta=delta1) 
  #
  # The index of locations in tc0
  grandstage <- rep(0, times=stageK1)
  for (i in 1:stageK1) {
    grandstage[i] <- utils::tail(which( tc0 <= tc1[i]), n=1)
  }
  # Debugging codes
  # tc0 <-c (3,6,9,12); tc1 <- c(3,4,7,9); stageK1 <- length(t1c); 
  # tc0 <-c (3,6,9,12); tc1 <- c(6,12,18,24); stageK1 <- length(t1c); 
  # print(grandstage)
  #
  # grandstage is jx #
  #
  finalgrandstage <- max(grandstage)
  # If the indices in grandstage are all the same, the decision on H1 doesn't depend on H0
  if (finalgrandstage == grandstage[1]) {
    pw1 <- pw1max
    return (pw1)
  } # End of if
  #
  subgrandstage <- grandstage[which(grandstage < finalgrandstage)]
  #
  subK0 <- length(subgrandstage)
  sqt0 <- sqrt(t0)
  sqt1 <- sqrt(t1)
  temp_result <- rep(0,times=subK0)
  for (i in 1:subK0) {
    gst0 <- subgrandstage[i] # grand stage
    gst1 <- i # grand stage
    meanV <- c(sqt0[1:gst0]*delta0, sqt1[gst1:stageK1]*delta1)
    corrM <- corrMatGenerator(tp=t0[1:gst0], ts=t1[gst1:stageK1], rhops=rho)
    lowerB <- c(rep(-Inf,times=gst0), cvec1[gst1], rep(-Inf,times=stageK1-gst1))
    upperB <- c(cvec0[1:gst0], Inf, cvec1[(gst1+1):stageK1])
    tempIntgl <- mvtnorm::pmvnorm(lowerB,upperB,meanV,corrM,algorithm=Miwa(steps=128))
    temp_result[i] <- tempIntgl[1]
  }
  pw1 <- pw1max - sum(temp_result)
  return (pw1)
}#
