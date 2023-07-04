# Rfun_sBoundsPh2
# This is different from sBndyPhrho or scdryBndyPhrho (2017-12-08)
# 2020-03-04/05
#
#' @name sBoundsPh2
#' @title Critical boundaries for the secondary endpoint in the partially hierarchical group sequential design
#' @description This function computes the critical boundaries for the secondary endpoint in the partially hierarchical group sequential design.
#' @param alpha a value of overall type I error rate
#' @param alpha0 a value of the type I error rate for testing H0
#' @param t0 a vector of information times for H0
#' @param t1 a vector of information times for H1
#' @param tc0 a vector of calendar times for H0
#' @param tc1 a vector of calendar times for H1
#' @param rho a value of the correlation between the test statistics for H0 and H1.
#' @param iuse0 a value of the type of error spending function for H0, ranging from 1 to 4
#' @param iuse1h a value of the type of error spending function for H1 first half, ranging from 1 to 4
#' @param iuse1t a value of the type of error spending function for H1 second half, ranging from 1 to 4
#' @param phi0 a value of the parameter of error spending function for H0
#' @param phi1h a value of the parameter of error spending function for H1 first half
#' @param phi1t a value of the parameter of error spending function for H1 second half
#' @return a list of two vectors: \code{bd} critical boundaries, \code{er} error spent until each stage
#' @export
#' @author Jiangtao Gou
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
#' sBoundsPh2(alpha, alpha0,
#'     t0, t1, tc0, tc1,
#'     rho, iuse0, iuse1h, iuse1t,
#'     phi0, phi1h, phi1t)
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
#' rho <- 0.0;
#' sBoundsPh2(alpha, alpha0,
#'     t0, t1, tc0, tc1,
#'     rho, iuse0, iuse1h=iuse1, iuse1t=iuse1,
#'     phi0, phi1h=phi1, phi1t=phi1)
#' @references
#'  Gou, J. (2023). Trigger strategy in repeated tests on multiple hypotheses. \emph{Statistics in Biopharmaceutical Research}, 15(1), 133-140.
#'  Gou, J. (2022). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. \emph{Biometrical Journal}, 64(2), 301-311.
#
#
sBoundsPh2 <- function(alpha, alpha0, t0, t1, tc0=t0, tc1=t1, rho=0, iuse0=1, iuse1h=1, iuse1t=1, phi0=rep(1,length(alpha)), phi1h=rep(1,length(alpha)), phi1t=rep(1,length(alpha))) {
  # Critical boundaries of H0
  cvecList0 <- gbounds(t=t0, iuse=iuse0, alpha=alpha0, phi=phi0)
  cvec0 <- cvecList0$bd
  #
  stageK0 <- length(t0)
  stageK1 <- length(t1)
  # Testing H0 lasts longer than testing H1
  if (is.na(which(tc1 >= tc0[stageK0])[1])) {
    idx_t1_tail <- stageK1
    stageK1_tail <- 1
    t1_tail <- 1
  } else {
    # Testing H1 lasts longer than testing H0
    idx_t1_tail <- which(tc1 >= tc0[stageK0]) # Pick the stages of H1 which are not skipped
    stageK1_tail <- length(idx_t1_tail)
    t1_tail <- t1[idx_t1_tail]
  } # End of if
  #
  # Part I Let the type I error of testing H0 intersect H1 less than alpha
  # alpha1t: marginal type I error for tail critical boundaries for testing H1
  targetDiff <- function(alpha1t, alpha, cvec0, t0, t1, tc0, tc1, rho, t1_tail, iuse1t, phi1t) {
    cvecList1t <- gbounds(t=t1_tail, iuse=iuse1t, alpha = alpha1t, phi=phi1t)
    cvec1t <- cvecList1t$bd
    # Head elembent of cvec1t_full are arbitrary
    cvec1t_full <- c(rep(0, times=length(t1)-length(t1_tail)), cvec1t)
    alpha01 <- sErrRphInt2(cvec0=cvec0, cvec1=cvec1t_full, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho)
    return(alpha-alpha01)
  } # End of Function targetDiff
  # lower bound and upper bound can be alpha-alpha0 and alpha
  alpha1tResult <- stats::uniroot(targetDiff, lower=(alpha-alpha0)/2, upper=2*alpha, tol=2.5e-16, alpha=alpha, cvec0=cvec0, t0=t0, t1=t1, tc0=tc0, tc1=tc1, rho=rho, t1_tail=t1_tail, iuse1t=iuse1t, phi1t=phi1t)
  alpha1t <- alpha1tResult$root # marginal type I error for tail of testing H1
  #
  # Part II Let the marginal type I error of cvec1 is less than alpha
  # In the extreme case, Delta0 is super larger, and H0 will be rejected immediately at the very beginning. Because of this, we need to control the marginial type I error of testing H1 at level alpha
  cvecList1t <- gbounds(t=t1_tail, iuse=iuse1t, alpha = alpha1t, phi=phi1t) # We have got alpha1t in Part I
  cvec1t <- cvecList1t$bd
  #
  targetDiff1 <- function(alpha1h, alpha, t1_head, t1_tail, iuse1h, phi1h, cvec1t){
    cvecList1h <- gbounds(t=t1_head, iuse=iuse1h, alpha = alpha1h, phi=phi1h)
    cvec1h <- cvecList1h$bd
    cvec1 <- c(cvec1h, cvec1t)
    t1 <- c(t1_head, t1_tail)
    alphas <- boundary2alpha(cvec=cvec1, t=t1)
    alpha_calc <- alphas[length(alphas)]
    return(alpha - alpha_calc)
  } # End of Function targetDiff1
  #
  t1_head <- t1[-idx_t1_tail]
  #
  alpha1hResult <- stats::uniroot(targetDiff1, lower=(alpha-alpha1t)/2, upper=2*alpha, tol=2.5e-16, alpha=alpha, t1_head=t1_head, t1_tail=t1_tail, iuse1h=iuse1h, phi1h=phi1h, cvec1t=cvec1t)
  alpha1h <- alpha1hResult$root
  #
  cvecList1h <- gbounds(t=t1_head, iuse=iuse1h, alpha = alpha1h, phi=phi1h) # We have got alpha1h above
  cvec1h <- cvecList1h$bd
  #
  cvec1 <- c(cvec1h, cvec1t)
  alphas1 <-boundary2alpha(cvec=cvec1, t=t1)
  #
  result <- list(bd=cvec1, er=alphas1)
  return(result)
}
