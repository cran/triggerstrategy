# Rfun_solveAlphaXsampleSize
# t0=1, t1=1 case: 2020-03-18 08:50
# trigger case: 2020-03-20 08:55
# other cases: 2020-03-20 11:55
# debug for other cases: 2020-03-21 12:50
#
#' @name solveAlphaXsampleSize
#' @title Sample size calculation
#' @description This function computes the sample size and the error rate pre-assigned to the primary endpoint using methods of \code{trigger}, \code{holm}, \code{maurer-bretz}, \code{bonferroni}.
#' @param alpha a number of overall type I error rate
#' @param beta0 a number of type II error rate for H0
#' @param beta1 a number of type II error rate for H1
#' @param effsz0 a number of the effect size of testing H0
#' @param effsz1 a number of the effect size of testing H1
#' @param szratio a number of the ratio of sample size of testing H0 to that of testing H1
#' @param t0 a vector of information times for H0
#' @param t1 a vector of information times for H1
#' @param tc0 a vector of calendar times for H0
#' @param tc1 a vector of calendar times for H1
#' @param rho a value of correlation coefficient between H0 and H1
#' @param iuse0 an integer shows the type of group sequential boundaries used for the primary endpoint
#' @param iuse1 an integer shows the type of group sequential boundaries used for the secondary endpoint
#' @param phi0 a parameter for the power family or the HSD gamma family for the primary endpoint
#' @param phi1 a parameter for the power family or the HSD gamma family for the secondary endpoint
#' @param usingRhoForBoundary an indicator whether using the informaiton of rho to calculate the boundary, default is FALSE (not using)
#' @param method a text of method, including \code{trigger}, \code{holm}, \code{maurer-bretz}, \code{bonferroni}
#' @param myinit a vector of two starting points for alpha0 and sample size.
#' @return a list of two values, \code{alpha0} and \code{groupsize}
#' @export
#' @import nleqslv
#' @import stats
#' @examples 
#' # Single Stage Example
#' alpha <- 0.025
#' effsz0 <- 0.4
#' effsz1 <- 0.30
#' szratio <- 1
#' beta0 <- 0.10
#' beta1 <- 0.20
#' solveAlphaXsampleSize(alpha, beta0, beta1, 
#'     effsz0, effsz1, szratio)
#' # Multi-stage example
#' alpha=0.025
#' beta0=0.10
#' beta1=0.20
#' effsz0=0.33
#' effsz1=0.30
#' szratio=1
#' t0=c(0.5,0.9,1)
#' t1=c(0.6,1)
#' tc0=c(1,2)
#' tc1=c(1,2,3)
#' rho=0
#' iuse0=1
#' iuse1=2
#' phi0=-4
#' phi1=1
#' usingRhoForBoundary=FALSE
#' myinit=c(300,alpha/2)
#' myinit=c(200,alpha/10)
#' method="trigger"
#' method="bonferroni"
#' method="holm"
#' method="maurer-bretz"
#' solveAlphaXsampleSize(alpha=alpha, 
#'     beta0=beta0, beta1=beta1, 
#'     effsz0=effsz0, effsz1=effsz1, 
#'     szratio=szratio, 
#'     t0=t0, t1=t1, tc0=tc0, tc1=tc1, 
#'     rho=rho, iuse0=iuse0, iuse1=iuse1, 
#'     phi0=phi0, phi1=phi1, 
#'     usingRhoForBoundary=usingRhoForBoundary, 
#'     method=method, 
#'     myinit=myinit) 
#' @references
#' Gou, J. (2021). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. Technical report. 
#
solveAlphaXsampleSize <- function(alpha, beta0, beta1, effsz0, effsz1, szratio=1, t0=1, t1=1, tc0=t0, tc1=t1, rho=0, iuse0=1, iuse1=1, phi0=rep(1,length(alpha)), phi1=rep(1,length(alpha)), usingRhoForBoundary=FALSE, method="trigger", myinit) {
  #
  alpha0 <- 0
  groupsize <- rep(0, times=2)
  #
  #<https://www.quora.com/What-is-the-difference-between-and-and-between-and-in-R>
  if (length(t0) == 1 && length(t1) == 1) {
    target1 <- function(alpha0, alpha, beta0, beta1, effsz0, effsz1) {
      lhs <- (qnorm(1-alpha0) + qnorm(1-beta0)) / (qnorm(1-alpha+alpha0) + qnorm(1-beta1))
      rhs <- effsz0/effsz1
      return(lhs-rhs)
    } # End of target
    result <- stats::uniroot(target1, lower=1e-9, upper=alpha-1e-9, tol=2.5e-16, alpha=alpha, beta0=beta0, beta1=beta1, effsz0=effsz0, effsz1=effsz1)
    alpha0 <- result$root
    groupsize[1] <- ((qnorm(1-alpha0) + qnorm(1-beta0))/effsz0)^2*(1+szratio)/szratio
    groupsize[2] <- szratio*groupsize[1]
    #
    result <- list(alpha0=alpha0, groupsize=groupsize)
    return(result)
  } # End of if
  #
  #
  target2 <- function(x, alpha, beta0, beta1, effsz0, effsz1, szratio, t0, t1, tc0, tc1, rho, iuse0, iuse1, phi0, phi1, usingRhoForBoundary) {
    groupsize <- x[1] 
    alpha0 <- x[2]
    #
    if (method == "trigger") {
      pspwr <- psPwRtrigger(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1, usingRhoForBoundary=usingRhoForBoundary, groupsize=groupsize, szratio=szratio, effsz0=effsz0, effsz1=effsz1) 
    } else if (method == "maurer-bretz") {
      pspwr <- psPwRbhmb(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1, usingRhoForBoundary=usingRhoForBoundary, groupsize=groupsize, szratio=szratio, effsz0=effsz0, effsz1=effsz1, method=method) 
    } else if (method == "bonferroni") {
      pspwr <- psPwRbhmb(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1, usingRhoForBoundary=usingRhoForBoundary, groupsize=groupsize, szratio=szratio, effsz0=effsz0, effsz1=effsz1, method=method) 
    } else if (method == "holm") {
      pspwr <- psPwRbhmb(alpha=alpha, alpha0=alpha0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1, usingRhoForBoundary=usingRhoForBoundary, groupsize=groupsize, szratio=szratio, effsz0=effsz0, effsz1=effsz1, method=method) 
    } else {
      stop("Methods include: bonferroni, holm, maurer-bretz, trigger.")
    } # End of if method
    y <- numeric(2)
    y[1] <- pspwr[1] - 1 + beta0
    y[2] <- pspwr[2] - 1 + beta1
    return(y)
  } # End of function target
  #
  xstart <- matrix(myinit,nrow=1,byrow=TRUE)
  ans <- nleqslv::searchZeros(xstart, target2, global="dbldog", alpha=alpha, beta0=beta0, beta1=beta1, effsz0=effsz0, effsz1=effsz1, szratio=szratio, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, iuse0=iuse0, iuse1=iuse1, phi0=phi0, phi1=phi1, usingRhoForBoundary=usingRhoForBoundary)
  #
  groupsize[1] <- ans$x[1]
  alpha0 <- ans$x[2]
  groupsize[2] <- szratio*groupsize[1]
  #
  result <- list(alpha0=alpha0, groupsize=groupsize)
  return(result)
}#


