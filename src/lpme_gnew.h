#ifndef _lpme_gnew_h
#define _lpme_gnew_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "lpme_common.h"

// estimate g function by ignoring measurement error using local polynomial
RcppExport SEXP fitlocpoly( SEXP x_, SEXP beta_, SEXP Kinput_, SEXP X_, SEXP Y_, SEXP h_);

// estimate g function when error is laplace
RcppExport SEXP fitnewLap( SEXP x_, SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP sigU_, SEXP h_);

// estimate g function when error is Gaussian
RcppExport SEXP fitnewGau( SEXP x_, SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP sigU_, SEXP h_);

// SIMEX bandwidth selection when error is laplace
RcppExport SEXP SIMEXnewLap( SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, SEXP pW_, SEXP pWs_);

// SIMEX bandwidth selection when error is Gaussian
RcppExport SEXP SIMEXnewGau( SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, SEXP pW_, SEXP pWs_);
#endif
