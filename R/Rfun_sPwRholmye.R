# Rfun_sPwRholmye
# 2020-03-06/07
#
#' @name sPwRholmye
#' @title Power of testing the secondary hypothesis using Holm-Ye
#' @description This function computes the power of testing the secondary hypothesis using Holm-Ye
#' @param alpha a number shows the overall error rate
#' @param alpha0 a number shows the error rate assigned to the primary endpoint initially
#' @param t0 a vector shows the information times of the primary endpoint
#' @param t1 a vector shows the information times of the secondary endpoint
#' @param tc0 a vector shows the calendar times of the primary endpoint
#' @param tc1 a vector shows the calendar times of the secondary endpoint
#' @param rho a number shows the correlation between the primary and secondary endpoints
#' @param iuse0 an integer shows the type of group sequential boundaries used for the primary endpoint
#' @param iuse1 an integer shows the type of group sequential boundaries used for the secondary endpoint
#' @param phi0 a parameter for the power family or the HSD gamma family for the primary endpoint
#' @param phi1 a parameter for the power family or the HSD gamma family for the secondary endpoint
#' @param delta0 a value of delta for hypothesis H0
#' @param delta1 a value of delta for hypothesis H1
#' @return a number shows the statistical power of rejecting H1
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
#' sPwRholmye(alpha=alpha, alpha0=alpha0, 
#'     t0=t0, t1=t1, tc0=tc0, tc1=tc1, 
#'     delta0=delta0, delta1=delta1, 
#'     rho=rho, iuse0=iuse0, iuse1=iuse1, 
#'     phi0=phi0, phi1=phi1) 
#
sPwRholmye <- function(alpha, alpha0, t0, t1, tc0=t0, tc1=t1, delta0, delta1, rho=0, iuse0=1, iuse1=1, phi0=rep(1,length(alpha)), phi1=rep(1,length(alpha))) {
  #
  stageK0 <- length(t0)
  stageK1 <- length(t1)
  #
  alpha1 <- alpha - alpha0
  #
  jx <- jxCalendarTime(tc0=tc0, tc1=tc1)
  jy <- jyCalendarTime(jx=jx)
  # print(jx); print(jy)
  #
  cvecList00 <- gbounds(t=t0, iuse=iuse0, alpha=alpha0, phi=phi0)
  cvec00 <- cvecList00$bd
  #
  cvecList11 <- gbounds(t=t1, iuse=iuse1, alpha=alpha1, phi=phi1)
  cvec11 <- cvecList11$bd
  cvecList1a <- gbounds(t=t1, iuse=iuse1, alpha=alpha, phi=phi1)
  cvec1a <- cvecList1a$bd
  #
  result_temp <- rep(0, times=length(jy))
  result_temp_id <- 0
  #
  for (k in jy) {
    if (k == jy[1]) {
      #
      result_temp_id <-  result_temp_id + 1
      #
      idy <- k:stageK1 #print(idy)
      meanV <- sqrt(t1[idy]) * delta1 #print(meanV)
      sigmaM <- corrMatGenerator(tp=vector(mode="numeric",length=0), ts=t1[idy], rhops=rho) #print(sigmaM)
      lowerB <- rep(-Inf,times=length(idy))
      upperB <- cvec11[idy]
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part1 <- resultIntgl[1] # print(Part1)
      #
      upperB <- cvec1a[idy]
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part2 <- resultIntgl[1] # print(Part1)
      #
      idx <- 1:jx[k]
      idy <- k:stageK1 
      meanV <- c(sqrt(t0[idx]) * delta0, sqrt(t1[idy]) * delta1) 
      sigmaM <- corrMatGenerator(tp=t0[idx], ts=t1[idy], rhops=rho) 
      lowerB <- rep(-Inf,times=length(idx)+length(idy))
      upperB <- c(cvec00[idx], cvec11[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part3 <- resultIntgl[1] 
      #
      upperB <- c(cvec00[idx], cvec1a[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part4 <- resultIntgl[1] 
      #
      # print(c(Part1, Part2, Part3, Part4))
      result_temp[result_temp_id] <- Part1 - Part2 - Part3 + Part4
    } else {
      #
      result_temp_id <-  result_temp_id + 1
      #
      idx <- 1:jx[k-1]
      idy <- k:stageK1 
      meanV <- c(sqrt(t0[idx]) * delta0, sqrt(t1[idy]) * delta1) 
      sigmaM <- corrMatGenerator(tp=t0[idx], ts=t1[idy], rhops=rho) 
      lowerB <- rep(-Inf,times=length(idx)+length(idy))
      upperB <- c(cvec00[idx], cvec11[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part1 <- resultIntgl[1] # print(Part1)
      #
      upperB <- c(cvec00[idx], cvec1a[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part2 <- resultIntgl[1] # print(Part1)
      #
      idx <- 1:jx[k]
      idy <- k:stageK1 
      meanV <- c(sqrt(t0[idx]) * delta0, sqrt(t1[idy]) * delta1) 
      sigmaM <- corrMatGenerator(tp=t0[idx], ts=t1[idy], rhops=rho) 
      lowerB <- rep(-Inf,times=length(idx)+length(idy))
      upperB <- c(cvec00[idx], cvec11[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part3 <- resultIntgl[1] 
      #
      upperB <- c(cvec00[idx], cvec1a[idy])
      resultIntgl <- mvtnorm::pmvnorm(lower=lowerB, upper=upperB, mean=meanV, sigma=sigmaM, algorithm=Miwa(steps=128))
      Part4 <- resultIntgl[1] 
      #
      result_temp[result_temp_id] <- Part1 - Part2 - Part3 + Part4
    }
  }
  pwr1Additional <- sum(result_temp)
  #
  pwr1Bonf <- sPwRnaiveBonf(alpha=alpha, alpha0=alpha0, t1=t1, delta1=delta1, iuse1=iuse1, phi1=phi1)
  #
  pwr1 <- pwr1Bonf + pwr1Additional
  #
  return(pwr1)
}