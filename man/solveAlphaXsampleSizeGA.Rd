% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rfun_solveAlphaXsampleSizeGA.R
\name{solveAlphaXsampleSizeGA}
\alias{solveAlphaXsampleSizeGA}
\title{Sample size calculation using Genetic Algorithms}
\usage{
solveAlphaXsampleSizeGA(
  alpha,
  beta0,
  beta1,
  effsz0,
  effsz1,
  szratio = 1,
  t0 = 1,
  t1 = 1,
  tc0 = t0,
  tc1 = t1,
  rho = 0,
  iuse0 = 1,
  iuse1 = 1,
  phi0 = rep(1, length(alpha)),
  phi1 = rep(1, length(alpha)),
  usingRhoForBoundary = FALSE,
  method = "trigger",
  lower = c(1, 1e-04),
  upper = c(10000, alpha - 1e-04),
  maxiter = 20,
  run = 200,
  seed = 1949
)
}
\arguments{
\item{alpha}{a number of overall type I error rate}

\item{beta0}{a number of type II error rate for H0}

\item{beta1}{a number of type II error rate for H1}

\item{effsz0}{a number of the effect size of testing H0}

\item{effsz1}{a number of the effect size of testing H1}

\item{szratio}{a number of the ratio of sample size of testing H0 to that of testing H1}

\item{t0}{a vector of information times for H0}

\item{t1}{a vector of information times for H1}

\item{tc0}{a vector of calendar times for H0}

\item{tc1}{a vector of calendar times for H1}

\item{rho}{a value of correlation coefficient between H0 and H1}

\item{iuse0}{an integer shows the type of group sequential boundaries used for the primary endpoint}

\item{iuse1}{an integer shows the type of group sequential boundaries used for the secondary endpoint}

\item{phi0}{a parameter for the power family or the HSD gamma family for the primary endpoint}

\item{phi1}{a parameter for the power family or the HSD gamma family for the secondary endpoint}

\item{usingRhoForBoundary}{an indicator whether using the informaiton of rho to calculate the boundary, default is FALSE (not using)}

\item{method}{a text of method, including \code{trigger}, \code{holm}, \code{maurer-bretz}, \code{bonferroni}}

\item{lower}{a vector of two lower limits for alpha0 and sample size}

\item{upper}{a vector of two upper limits for alpha0 and sample size.}

\item{maxiter}{a number of maximum number of iterations}

\item{run}{a number of maximum number of consecutive generations without any improvement in the best fitness value before the GA is stopped}

\item{seed}{a number of seed of random number generator}
}
\value{
a list of two values, \code{alpha0} and \code{groupsize}
}
\description{
This function computes the sample size and the error rate pre-assigned to the primary endpoint using methods of \code{trigger}, \code{holm}, \code{maurer-bretz}, \code{bonferroni}, with Genetic Algorithms.
}
\details{
R package \code{GA} is used for Genetic Algorithms.
}
\examples{
alpha=0.025
beta0=0.10
beta1=0.20
effsz0=0.33
effsz1=0.30
szratio=1
t0=c(0.5,0.9,1)
t1=c(0.6,1)
tc0=c(1,2)
tc1=c(1,2,3)
rho=0
iuse0=1
iuse1=2
phi0=-4
phi1=1
usingRhoForBoundary=FALSE
method="trigger"
method="bonferroni"
method="holm"
method="maurer-bretz"
lower = c(180,0.005)
upper = c(240, alpha-0.005)
maxiter = 1 # Increase this number for more precise results
run = 1 # Increase this number for more precise results
seed = 123
result <- solveAlphaXsampleSizeGA(alpha=alpha,
    beta0=beta0, beta1=beta1,
    effsz0=effsz0, effsz1=effsz1,
    szratio=szratio,
    t0=t0, t1=t1, tc0=tc0, tc1=tc1,
    rho=rho, iuse0=iuse0, iuse1=iuse1,
    phi0=phi0, phi1=phi1,
    usingRhoForBoundary=usingRhoForBoundary,
    method=method,
    lower = lower, upper = upper,
    maxiter = maxiter,
    run = run,
    seed = seed)
print(result)
}
\references{
Gou, J. (2023). Trigger strategy in repeated tests on multiple hypotheses. \emph{Statistics in Biopharmaceutical Research}, 15(1), 133-140.
Gou, J. (2022). Sample size optimization and initial allocation of the significance levels in group sequential trials with multiple endpoints. \emph{Biometrical Journal}, 64(2), 301-311.
}
