# Rfun_psPwRbhmb
# 2020-03-20 11:30
#
#' @name psPwRbhmb
#' @title Powers of testing the primary and secondary hypotheses
#' @description This function computes the powers of testing the primary and secondary hypotheses using the \code{holm}, \code{maurer-bretz}, \code{bonferroni} methods. 
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
#' @param usingRhoForBoundary an indicator whether using the informaiton of rho to calculate the boundary, default is FALSE (not using)
#' @param groupsize a value of sample size in group 1
#' @param szratio a value of the sample size ratio, n2/n1
#' @param effsz0 a value of effect size for hypothesis H0
#' @param effsz1 a value of effect size for hypothesis H1
#' @param method a text of method, including \code{holm}, \code{maurer-bretz}, \code{bonferroni}
#' @return a vector of two values of the probability that H0 is rejected,  the probability that H1 is rejected, the power, using \code{holm}, \code{maurer-bretz}, or \code{bonferroni}. 
#' @export
#' @import stats
#' @author Jiangtao Gou
#' @details The methods include \code{holm}, \code{maurer-bretz}, and \code{bonferroni}. Users can decide whether the correlation information is used or not.
#' @examples 
#' alpha <-  0.025
#' alpha0 <- 0.0136
#' iuse0 <- 1
#' iuse1 <- 2
#' phi0 <- -4
#' phi1 <- 1
#' tc0 <- c(1,2)
#' tc1 <- c(1,2,3)
#' t0 <- c(0.6,1)
#' t1 <- c(0.5,0.9,1)
#' rho <- 0
#' effsz0 <- 0.33 
#' effsz1 <- 0.30
#' groupsize=226
#' szratio=1
#' method="bonferroni"
#' method = "holm"
#' method="maurer-bretz"
#' psPwRbhmb(alpha=alpha, alpha0=alpha0, 
#'     t0=t0, t1=t1, tc0=tc0, tc1=tc1, 
#'     rho=rho, iuse0=1, iuse1=iuse1, 
#'     phi0=phi0, phi1=phi1, 
#'     usingRhoForBoundary=usingRhoForBoundary, 
#'     groupsize=groupsize, szratio=szratio, 
#'     effsz0=effsz0, effsz1=effsz1, 
#'     method=method) 
#' @references
#' Gou, J. (2021). Trigger strategy in repeated tests on multiple hypotheses. Technical report. 
#
psPwRbhmb <- function(alpha, alpha0, t0, t1, tc0=t0, tc1=t1, rho=0, iuse0=1, iuse1=1, phi0=rep(1,length(alpha)), phi1=rep(1,length(alpha)), usingRhoForBoundary=FALSE, groupsize, szratio=1, effsz0, effsz1, method="bonferroni") {
  #
  pspwr <- rep(0, times=2)
  #
  n1 <- groupsize
  n2 <- szratio*groupsize
  delta0 <- effsz0*sqrt(n1*n2/(n1+n2))
  delta1 <- effsz1*sqrt(n1*n2/(n1+n2))
  #
  if (method == "bonferroni") {
    pspwr[1] <- sPwRnaiveBonf(alpha=alpha, alpha0=alpha-alpha0, t1=t0, delta1=delta0, iuse1=iuse0, phi1=phi0)
    pspwr[2] <- sPwRnaiveBonf(alpha=alpha, alpha0=alpha0, t1=t1, delta1=delta1, iuse1=iuse1, phi1=phi1)
  } else if (method == "holm") {
    pspwr[1] <- sPwRholm(alpha=alpha, alpha0=alpha-alpha0, t0=t1, t1=t0, delta0=delta1, delta1=delta0, rho=rho, iuse0=iuse1, iuse1=iuse0, phi0=phi1, phi1=phi0)
    pspwr[2] <- sPwRholm(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, delta0=delta0, delta1=delta1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1)
  } else if (method == "maurer-bretz") {
    pspwr[1] <- sPwRholmye(alpha=alpha, alpha0=alpha-alpha0, t0=t1, t1=t0, tc0=tc1, tc1=tc0, delta0=delta1, delta1=delta0, rho=rho, iuse0=iuse1, iuse1=iuse0, phi0=phi1, phi1=phi0) 
    pspwr[2] <- sPwRholmye(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, delta0=delta0, delta1=delta1, rho=rho, iuse0=1, iuse1=iuse1, phi0=phi0, phi1=phi1) 
  } else {
    stop("Methods include: bonferroni, holm, maurer-bretz.")
  } # End of if method
  #
  return (pspwr)
}