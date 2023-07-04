# Rfun_sPwRtrigger
# 2020-03-17
#
#' @name sPwRtrigger
#' @title Power of testing the secondary hypothesis using Trigger strategy
#' @description This function computes the power of testing the secondary hypothesis using  Trigger strategy
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
#' @param usingRhoForBoundary an indicator whether using the informaiton of rho to calculate the boundary, default is FALSE (not using)
#' @return a value of the probability that H1 is rejected, the power, using the trigger strategy
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
#' sPwRtrigger(alpha=alpha, alpha0=alpha0,
#'     t0=t0, t1=t1, tc0=tc0, tc1=tc1,
#'     delta0=delta0, delta1=delta1,
#'     rho=rho, iuse0=1, iuse1=iuse1,
#'     phi0=phi0, phi1=phi1,
#'     usingRhoForBoundary=FALSE)
#' @references
#'  Gou, J. (2023). Trigger strategy in repeated tests on multiple hypotheses. \emph{Statistics in Biopharmaceutical Research}, 15(1), 133-140.
#'  Gou, J. (2022). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. \emph{Biometrical Journal}, 64(2), 301-311.
#
sPwRtrigger <- function(alpha, alpha0, t0, t1, tc0=t0, tc1=t1, delta0, delta1, rho=0, iuse0=1, iuse1=1, phi0=rep(1,length(alpha)), phi1=rep(1,length(alpha)), usingRhoForBoundary=FALSE) {
  cvecList0 <- gbounds(t=t0, iuse=iuse0, alpha=alpha0, phi=phi0)
  cvec0 <- cvecList0$bd
  if (!usingRhoForBoundary) {
  # Boundaries are good for all non-negative correlations
  cvecList1 <- sBoundsPh2(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=0, iuse0=iuse0, iuse1h=iuse1, iuse1t=iuse1, phi0=phi0, phi1h=phi1, phi1t=phi1)
  } else {
    # Boundaries are good for a specified correlation
    cvecList1 <- sBoundsPh2(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1h=iuse1, iuse1t=iuse1, phi0=phi0, phi1h=phi1, phi1t=phi1)
  }
  cvec1 <- cvecList1$bd
  pwr1 <- sPwRph2(cvec0=cvec0, cvec1=cvec1, delta0=delta0, delta1=delta1, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho)
  return(pwr1)
}
