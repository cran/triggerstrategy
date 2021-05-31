# Rfun_boundary2alpha
#
#' @name boundary2alpha
#' @title Convert normal critical boundaries to cumulative alpha levels. 
#' @description This function converts normal critical boundaries to cumulative alpha levels.
#' @param cvec a vector of critical boundaries
#' @param t a vector of information times
#' @return alphas, a vector of cumulative errors
#' @export
#' @author Jiangtao Gou
#' @author Fengqing Zhang
#' @examples
#' t <- c(0.5,0.8,1)
#' iuse <- 4
#' result <- gbounds(t=t, iuse=iuse)
#' print(result)
#' boundary2alpha(cvec=result$bd, t=t)
#
boundary2alpha <- function(cvec, t) {
  K <-length(cvec)
  alphas <- rep(0, times=K)
  for (i in 1:K) {
    alphas[i] <- marginalPwR(cvec=cvec[1:i], t=t[1:i], delta=0)
  }
  return(alphas)
}