#ifndef _lpme_common_h
#define _lpme_common_h

#include <RcppArmadillo.h>

#define ESMALL 1e-10  /* small number */
#define ELARGE 1e+10 /* large number */
#define sqrt2pi 2.5066282746310002416 /* sqrt(2*pi)*/
typedef Rcpp::NumericMatrix::iterator mat_iterator;
using namespace Rcpp;

// Fourier transform for error U
NumericVector FuLap(NumericVector t, double sigU);
NumericVector FuGau(NumericVector t, double sigU);
NumericVector FuLapinv(NumericVector t, double sigU);
NumericVector FuGauinv(NumericVector t, double sigU);
// Fourier transform for Kernel K
NumericVector FK(NumericVector t);

// first derivative for Fourier transform of Kernel K
NumericVector FK1(NumericVector t);

// second derivative for Fourier transform of Kernel K
NumericVector FK2(NumericVector t);

// Fourier transform for Kernel K
NumericVector FK_sec_order(NumericVector t);
// first derivative for Fourier transform of Kernel K
NumericVector FK1_sec_order(NumericVector t);
// second derivative for Fourier transform of Kernel K
NumericVector FK2_sec_order(NumericVector t);

// function to generate subvectors a[-(ind1:ind2)] and b[-(ind1:ind2)] and save them to w and y respectively
// void allows more than two values can be returned. For example, w and y are returned here.
void subvecij(const NumericVector& a, const NumericVector& b, int ind1, int ind2, NumericVector& w, NumericVector& y);

// function to estimate ghat of x using JASA when U is Laplace
void gjasaLap(NumericVector& res, const NumericVector& x, const NumericVector& t, double dt, const NumericVector& W, 
      const NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using JASA when U is Gaussian
void gjasaGau(NumericVector& res, const NumericVector& x, const NumericVector& t, double dt, const NumericVector& W, 
      const NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using OURS
void gnewLap(NumericVector& res, const NumericVector& x, const NumericVector& input, const NumericVector& output, 
        double beta, double beta2, const NumericVector& mconst, const NumericVector& Kinput, 
        const NumericVector& W, const NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using OURS when error is Gaussian
void gnewGau(NumericVector& ghatofx, const NumericVector& x, const NumericVector& input, const NumericVector& output, 
        double beta, double beta2, const NumericVector& mconst, const NumericVector& Kinput, 
        const NumericVector& W, const NumericVector& Y, double sigU, double h);

#endif
