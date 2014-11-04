#ifndef _lpme_gjasa_h
#define _lpme_gjasa_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "lpme_common.h"

// estimate g function JASA when error is laplace
RcppExport SEXP fitjasaLap( SEXP x_, SEXP h_, SEXP W_, SEXP Y_, SEXP sigU_, SEXP dt_, SEXP t_);

// estimate g function JASA when error is Gaussian
RcppExport SEXP fitjasaGau( SEXP x_, SEXP h_, SEXP W_, SEXP Y_, SEXP sigU_, SEXP dt_, SEXP t_);

// SIMEX bandwidth selection when error is laplace
RcppExport SEXP SIMEXjasaLap(SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, 
                  SEXP pW_, SEXP pWs_, SEXP dt_, SEXP t_);

// SIMEX bandwidth selection when error is Gaussian
RcppExport SEXP SIMEXjasaGau(SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, 
                  SEXP pW_, SEXP pWs_, SEXP dt_, SEXP t_);

#endif
