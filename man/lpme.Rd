\name{lpme}
\alias{lpme}
\alias{lpme.default}

\title{Local Polynomial Estimators for the Errors-in-Variables Problem}
\description{
This function provides both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2014+) local polynomial estimators for solving the errors-in-variables problem.
}
   
\usage{
lpme(Y, W, bw, method="HZ", sig=NULL, error="laplace", xgridmin=-2, xgridmax=2,
     data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{Y}{         an n by 1 response vector. }
   
   \item{W}{         an n by 1 predictor vector.} 
           
  \item{bw}{         bandwidth.} 
    
    \item{method}{   the method to be used; \code{method="HZ"} returns the estimator proposed by Huang and Zhou (2015); 
                     \code{method="DFC"} returns the estimator proposed by Delaigle, Fan, and Carroll (2009); \code{method="naive"}
                     returns the local polynomial estimator ignoring measurement error.}

    \item{sig}{      standard deviation of the measurement error.}   

    \item{error}{    the distribution assumed for the measurement error; \code{error="laplace"} is for Laplace distribution; 
                     \code{error="normal"} is for Gaussian distribution.} 
                    
    \item{xgridmin}{ the minimum value where the estimated responsed is evaluated at; default is -2.}
    
    \item{xgridmax}{ the maximum value where the estimated responsed is evaluated at; defualt is 2.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{lpme} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
This function provides both the DFC (Delaigle, Fan, and Carroll, 2009) and HZ (Huang and Zhou, 2014+) local polynomial estimators for solving the errors-in-variables problem.
}

\value{
The results include the grid points for predictor \code{xgrid} and corresponding fitted responses \code{yhat}.
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

## sample size:
n =100;
## Function gofx(x) to estimate
gofx  = function(x){ 1/4*x + x^2/4 }

## Generate data
sigma_e  = 0.5;
sigma_x = 1; X = rnorm(n, 0, sigma_x); 
## Sample Y
Y  = gofx(X) + rnorm(n, 0, sigma_e);
##------------------  method Based on X ---------------------------
#ghat_X= lpme(Y, X, 0.1, method="naive");

## reliability ratio
lambda=0.85;
sigma_u  = sqrt(1/lambda-1)*sigma_x;
print( sigma_x^2/(sigma_x^2 + sigma_u^2) );
W=X+rnorm(n,0,sigma_u);
#W=X+rlaplace(n,0,sigma_u/sqrt(2));
  
##------------------  method Based on W ---------------------------
#ghat_W=lpme(Y, W, 0.1, method="naive");

##------------------ JASA method -------------------------------------
h = 0.13; 
#ghat_JASA = lpme(Y, W, h, method="DFC", sig=sigma_u, error="laplace");

##------------------ Our method -------------------------
ghat_NEW = lpme(Y, W, h, method="HZ", sig=sigma_u, error="laplace");

## plots
x = ghat_NEW$xgrid;
plot(x, gofx(x), "l", main="Individual", lwd="2")
lines(ghat_NEW$xgrid, ghat_NEW$yhat, lty="dashed", col="2",lwd="3")
#lines(ghat_JASA$xgrid, ghat_JASA$yhat, lty="dotted", col="3",lwd="3")
#lines(ghat_X$xgrid, ghat_X$yhat, lty="dashed", col="4",lwd="2")
#lines(ghat_W$xgrid, ghat_W$yhat, lty="dashed", col="5",lwd="3")
}

