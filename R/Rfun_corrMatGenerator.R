# function corrMatGenerator 
# 2017-12-04
# 2020-03-07 updated
#
#' @name corrMatGenerator
#' @title Correlation matrix generator
#' @description This function generate the correlation matrix for group sequential trials with two endpoints.
#' @param tp an information fraction vector of Hp
#' @param ts an information fraction vector of Hs
#' @param rhops a number that shows the correlation coefficient between the test statistics of the primary endpoint and that of the secondary endpoint
#' @return the correlation matrix, Hp goes first, and Hs goes second. 
#' @export
#' @import stats
#' @author Jiangtao Gou
#' @examples
#' corrMatGenerator(tp=c(0.64,1),
#'     ts=c(0.25,0.49,1),
#'     rhops=1)
#
corrMatGenerator <- function (tp,ts,rhops) {
  # Either tp or ts can be NULL, 2020-03-07 18:25
  if (length(tp) == 0) {
    tp <- ts
    ts <- vector(mode="numeric",length=0)
  }
  Kp <- length(tp)
  Ks <- length(ts)
  K <- Kp + Ks
  lmd <- sqrt(tp)
  gmm <- sqrt(ts)
  corrMat <- matrix(rep(1,K*K),nrow=K,ncol=K,byrow=TRUE)
  for (i in 1:Kp) {
    for (j in 1:Kp){
      corrMat[i,j] <- min(lmd[i],lmd[j])/max(lmd[i],lmd[j])
    }
  }
  if (length(ts) == 0) {
    return(corrMat)
  }
  for (i in 1:Ks) {
    for (j in 1:Ks){
      corrMat[Kp+i,Kp+j] <- min(gmm[i],gmm[j])/max(gmm[i],gmm[j])
    }
  }
  for (i in 1:Kp) {
    for (j in 1:Ks) {
      corrMat[i,Kp+j] <- rhops * min(lmd[i],gmm[j])/max(lmd[i],gmm[j])
    }
  }
  for (i in 1:Ks) {
    for (j in 1:Kp) {
      corrMat[Kp+i,j] <- corrMat[j,Kp+i]
    }
  }
  return(corrMat)
}#
