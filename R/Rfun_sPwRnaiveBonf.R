# Rfun_sPwRnaiveBonf
# 2020-03-05
#
#' @name sPwRnaiveBonf
#' @title Power of testing the secondary hypothesis using Bonferroni
#' @description This function computes the power of testing the secondary hypothesis using Bonferroni
#' @param alpha a number shows the overall error rate
#' @param alpha0 a number shows the error rate assigned to the primary endpoint initially
#' @param t1 a vector shows the information times of the secondary endpoint
#' @param iuse1 an integer shows the type of group sequential boundaries used for the secondary endpoint
#' @param phi1 a parameter for the power family or the HSD gamma family for the secondary endpoint
#' @param delta1 a value of delta for hypothesis H1
#' @return a value of the probability that H1 is rejected, the power, using the naive Bonferroni method
#' @export
#' @examples
#' alpha <-  0.025
#' alpha0 <- 0.01
#' iuse0 <- 4
#' iuse1 <- 4
#' phi0 <- -4
#' phi1 <- -4
#' tc0 <- c(3,6,9,12,18)
#' tc1 <- c(6,12,18,36)
#' t0 <- (1:5)/5
#' t1 <- (1:4)/4
#' rho <- 0.5
#' delta0 <- 1
#' delta1 <- 3
#' sPwRnaiveBonf(alpha=alpha, 
#'     alpha0=alpha0, 
#'     t1=t1, 
#'     delta1=delta1, 
#'     iuse1=iuse1, 
#'     phi1=phi1)
#
sPwRnaiveBonf <- function(alpha, alpha0=alpha/2, t1, delta1, iuse1, phi1=rep(1,length(alpha))) {
  alpha1 <- alpha - alpha0
  cvecList1 <- gbounds(t=t1, iuse=iuse1, alpha=alpha1, phi=phi1)
  cvec1 <- cvecList1$bd
  pwr1 <- marginalPwR(cvec=cvec1,t=t1,delta=delta1)
  return(pwr1)
}
