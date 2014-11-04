\name{bwSIMEX}
\alias{bwSIMEX}
\alias{bwSIMEX.default}

\title{CV Bandwidth Selector Using SIMEX}
\description{
This function selects the bandwidth for both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2014+) estimators. 
}
   
\usage{
bwSIMEX(Y, W, method="HZ", sig, error="laplace", k_fold=2, B=10, 
        h1=NULL, h2=NULL,  length.h=10, Wdiff=NULL,
        data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
    
    \item{method}{   the method to be used; \code{method="HZ"} uses the estimator proposed by Huang and Zhou (2014+); 
                     \code{method="DFC"} uses the estimator proposed by Delaigle, Fan, and Carroll (2009).}

    \item{sig}{      standard deviation of the measurement error.}   

    \item{error}{    the distribution assumed for the measurement error; \code{error="laplace"} is for Laplace distribution; 
                     \code{error="normal"} is for Gaussian distribution.} 
                    
    \item{k_fold}{   gives fold of cross-validation to be used; default is 2.}
    
    \item{B}{        total number of cross-validation criteria to average over; defualt is 10.}
    
    \item{h1}{       bandwidth vector for the first level error contamination; default is \code{NULL}, and h1 is chosen automatically. See Huang and Zhou (2014+) for details. }
    
    \item{h2}{      bandwidth vector for the first level error contamination; defualt is \code{NULL}, and h2 is chosen automatically. See Huang and Zhou (2014+) for details.}
    
    \item{length.h}{ number of grid points for each of h1 and h2; default is 10.}
    
    \item{Wdiff}{    an n by 1 vector of \code{(W1-W2)/2}, where W1, W2 are two replicated measurements; default is \code{NULL}, which indicates that the errors are generated from the assumed error distribution, otherwise, the errors are generated from Wdiff with replacement.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{bwSIMEX} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
This function selects the bandwidth for both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2014+) estimators. 
}

\value{
  The results include the bandwidth \code{bw}.
}


\references{
Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice in errors-in-variables problems. \emph{Journal of the American Statistical Association}, 103, 280-287.

Delaigle, A., Fan, J., and Carroll, R. (2009). A design-adaptive local polynomial estimator for the errors-in-variables problem. \emph{Journal of the American Statistical Association} 104, 348-359.

Huang, X. and Zhou, H. (2014+). An alternative local polynomial estimator for the errors-in-variable problem. \emph{Submitted}.
}


\author{
Haiming Zhou <zhouh@email.sc.edu> and Xianzheng Huang <huang@stat.sc.edu>
}

\examples{
#############################################
## X - True covariates
## W - Observed covariates
## Y - individual response
rm(list=ls())
library(lpme)
## generate laplace
rlap=function (use.n, location = 0, scale = 1) 
{
location <- rep(location, length.out = use.n)
scale <- rep(scale, length.out = use.n)
rrrr <- runif(use.n)
location-sign(rrrr-0.5)*scale*(log(2)+ifelse(rrrr<0.5, log(rrrr), log1p(-rrrr)))
}

## sample size:
n =100;
## Function gofx(x) to estimate
gofx  = function(x){ 2*x*exp(-10*x^4/81) }

## Generate data
sigma_e  = 0.2;
sigma_x = 1; X = rnorm(n, 0, sigma_x); 
## Sample Y
Y  = gofx(X) + rnorm(n, 0, sigma_e);
## reliability ratio
lambda=0.85;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
print( sigma_x^2/(sigma_x^2 + sigma_u^2) );
#W=X+rnorm(n,0,sigma_u);
W=X+rlap(n,0,sigma_u/sqrt(2));
  
#### SIMEX
hwNEW = bwSIMEX(Y, W, method="HZ", sig=sigma_u, error="laplace", k_fold=2, 
                B=1, length.h=1)$bw
ghat_NEW = lpme(Y, W, hwNEW , method="HZ", sig=sigma_u, error="laplace");

## plots
x = ghat_NEW$xgrid;
plot(x, gofx(x), "l", main="Individual", lwd="2")
lines(ghat_NEW$xgrid, ghat_NEW$yhat, lty="dashed", col="2",lwd="3")
}

