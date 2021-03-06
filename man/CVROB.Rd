% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{CVROB}
\alias{CVROB}
\title{RoBoost-PLSR : Robust method for partial least squares regression}
\usage{
CVROB(X, Y, segm, alpha, beta, gamma, ncomp, niter, th)
}
\arguments{
\item{X}{Explanatory variables}

\item{Y}{Explained Variables}

\item{segm}{list of index}

\item{alpha}{parameters for X-residuals}

\item{beta}{parameters for Y-residuals}

\item{gamma}{parameters for leverage point}

\item{ncomp}{Number of latent variables}

\item{niter}{Number of maximal iterations}

\item{th}{convergence threshold}
}
\value{

}
\description{
RoBoost-PLSR : Robust method for partial least squares regression
}
\examples{


n <- 100
p <- 10
set.seed(1)
X <- matrix(rnorm(n * p, mean = 10), ncol = p)
Y <- 100 * rnorm(n)

library(dplyr)
n =nrow(X)
segm = segmkf(n = n,K = 5, nrep = 5)
ncomp = 5
th = 1-10^-12
niter = 20
alpha = Inf
beta = 4
gamma = Inf
test = CVROB(X = X, Y = Y ,segm = segm,alpha = alpha,beta = beta,gamma = gamma,ncomp = ncomp,niter = niter,th = th)
rmsecv = test \%>\% group_by(ncomp) \%>\% summarise(rmse.t(r,b,z = 0.15))
colnames(rmsecv)[2] = "r"
rmsecv[rmsecv$r==min(rmsecv$r),]
r2cv = test \%>\% group_by(ncomp) \%>\% summarise(r2.t(r,Yt,b,z = 0.15))
colnames(r2cv)[2] = "r"
r2cv[r2cv$r==max(r2cv$r),]
}
