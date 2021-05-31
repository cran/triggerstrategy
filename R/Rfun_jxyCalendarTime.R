# Rfun_jxyCalendarTime
# 2020-03-07
#
# Example 
#tc0 <- c(3,6,9,12); tc1 <- c(6,12,24)
#tc0 <- c(3,6,9,12); tc1 <- c(3,4,6,9)
#tc0 <- c(3,6,9,12); tc1 <- c(1,3,4,6,9)
#tc0 <- c(3,6,9,12,18); tc1 <- c(6,12,18,36);
#jx <- jxCalendarTime(tc0, tc1); print(jx)
#jy <- jyCalendarTime(jx); print(jy)
#
jxCalendarTime <- function(tc0, tc1) {
  jx <- rep(0, times=length(tc1))
  for (k in 1:length(tc1)) {
    idx <- which(tc0 <= tc1[k])
    if (length(idx) == 0) {
      jx[k] <- 0
    } else {
    jx[k] <- idx[length(idx)]
    }
  }
  return(jx)
}#

jyCalendarTime <- function(jx) {
  jx0 <- c(0, jx[1:(length(jx)-1)])
  jx <- which(jx > jx0)
  return(jx)
}#