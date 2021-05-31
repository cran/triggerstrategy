# Rfun_sPwRholm
# 2020-03-05
# Simple Holm, doesn't depend on calendar time
#
#' @name sPwRholm
#' @title Power of testing the secondary hypothesis using Holm
#' @description This function computes the power of testing the secondary hypothesis using Holm
#' @param alpha a number shows the overall error rate
#' @param alpha0 a number shows the error rate assigned to the primary endpoint initially
#' @param t0 a vector shows the information times of the primary endpoint
#' @param t1 a vector shows the information times of the secondary endpoint
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
#' t0 <- (1:5)/5
#' t1 <- (1:4)/4
#' rho <- 0.5
#' delta0 <- 1
#' delta1 <- 3
#' sPwRholm(alpha=alpha, alpha0=alpha0, 
#'     t0=t0, t1=t1, 
#'     delta0=delta0, delta1=delta1, 
#'     rho=rho, iuse0=iuse0, iuse1=iuse1, 
#'     phi0=phi0, phi1=phi1)
#
sPwRholm <- function(alpha, alpha0, t0, t1, delta0, delta1, rho=0, iuse0=1, iuse1=1, phi0=rep(1,length(alpha)), phi1=rep(1,length(alpha))){
  #
  stageK0 <- length(t0)
  stageK1 <- length(t1)
  #
  alpha1 <- alpha - alpha0
  #
  cvecList00 <- gbounds(t=t0, iuse=iuse0, alpha=alpha0, phi=phi0)
  cvec00 <- cvecList00$bd
  # cvecList0a <- gbounds(t=t0, iuse=iuse0, alpha=alpha, phi=phi0)
  # cvec0a <- cvecList0a$bd
  #
  cvecList11 <- gbounds(t=t1, iuse=iuse1, alpha=alpha1, phi=phi1)
  cvec11 <- cvecList11$bd
  cvecList1a <- gbounds(t=t1, iuse=iuse1, alpha=alpha, phi=phi1)
  cvec1a <- cvecList1a$bd
  # Part 1
  Part1 <- marginalPwR(cvec=cvec1a,t=t1,delta=delta1)
  # Part 2 & 3
  meanV <- c(sqrt(t0)*delta0, sqrt(t1)*delta1)
  corrM <- corrMatGenerator(tp=t0, ts=t1, rhops=rho)
  lowerB <- rep(-Inf,times=stageK0+stageK1)
  # Part 2
  upperB <- c(cvec00, cvec11)
  resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
  Part2 <- resultIntgl[1]
  # Part 3
  upperB <- c(cvec00, cvec1a)
  resultIntgl <- mvtnorm::pmvnorm(lowerB, upperB, meanV, corrM, algorithm=Miwa(steps=128))
  Part3 <- resultIntgl[1]
  #
  pwr1 <- Part1 - Part2 + Part3
  #print(c(Part1, Part2, Part3))
  #
  return(pwr1)
}

# Power of Graphical Approach
# \begin{align*}
# P_1 &= 1 - \Pr\left(  \{\cap_i \{Y < d_l\}\}\right) \\
# P_2 &= \Pr\left(\{\cap_i \{X < c\}\}^c \cap \{\cap_i \{Y < d_s\}\}^c \cap \{\cap_i \{Y < d_l\}\}\right) \\
# &= \Pr\left(\{\cap_i \{X < c\}\}^c \cap \{\cap_i \{Y < d_l\}\}\right) -  \Pr\left(\{\cap_i \{X < c\}\}^c \cap \{\cap_i \{Y < d_s\}\}\right) \\
# &= \Pr\left(  \{\cap_i \{Y < d_l\}\}\right) - \Pr\left(\{\cap_i \{X < c\}\} \cap \{\cap_i \{Y < d_l\}\}\right) \\
# &\qquad - \Pr\left(  \{\cap_i \{Y < d_s\}\}\right) + \Pr\left(\{\cap_i \{X < c\}\} \cap \{\cap_i \{Y < d_s\}\}\right)\\
# P &= P_1 + P_2 \\
# &= 1 - \Pr\left(\{\cap_i \{X < c\}\} \cap \{\cap_i \{Y < d_l\}\}\right) \\
# &\qquad - \Pr\left(  \{\cap_i \{Y < d_s\}\}\right) + \Pr\left(\{\cap_i \{X < c\}\} \cap \{\cap_i \{Y < d_s\}\}\right)
# \end{align*}