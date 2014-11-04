#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma ;
using namespace Rcpp ;
using namespace std;

// Fourier transform for error U
NumericVector FuLap(NumericVector t, double sigU){
  return(1.0/(1.0+Rcpp::pow((sigU*t),2)*0.5));
}
NumericVector FuGau(NumericVector t, double sigU){
  return( Rcpp::exp(-Rcpp::pow((sigU*t),2)*0.5) );
}

// Fourier transform for Kernel K
NumericVector FK(NumericVector t){
  return(Rcpp::pow((1.0-Rcpp::pow(t,2)), 8) );
}

// first derivative for Fourier transform of Kernel K
NumericVector FK1(NumericVector t){
  return(-16.0*t*Rcpp::pow((1.0-Rcpp::pow(t,2)), 7));
}

// second derivative for Fourier transform of Kernel K
NumericVector FK2(NumericVector t){
  return(-16.0*Rcpp::pow((1.0-Rcpp::pow(t,2)), 7) + 16.0*14.0*t*t*Rcpp::pow((1.0-Rcpp::pow(t,2)), 6));
}

// function to generate subvectors a[-(ind1:ind2)] and b[-(ind1:ind2)] and save them to w and y respectively
// void allows more than two values can be returned. For example, w and y are returned here.
void subvecij(const NumericVector& a, const NumericVector& b, int ind1, int ind2, NumericVector& w, NumericVector& y){
  int ny = y.size();
  for (int i=0; i<ny; ++i){
    if(i<ind1) {
      w[i] = a[i];
      y[i] = b[i];
    } else {
      w[i] = a[i+ind2-ind1+1];
      y[i] = b[i+ind2-ind1+1];
    }
  }
}

// function to estimate ghat of x using JASA when U is laplace
void gjasaLap(NumericVector& res, const NumericVector& x, const NumericVector& t, double dt, const NumericVector& W, 
      const NumericVector& Y, double sigU, double h){
  int nt = t.size();
  int n = W.size();
  int nx = x.size();
  NumericVector rehatFW(nt); 
  NumericVector imhatFW(nt);
  NumericVector reYhatFW(nt);
  NumericVector imYhatFW(nt);
  NumericVector FKt = FK(t);
  NumericVector FKt1= FK1(t);
  NumericVector FKt2= FK2(t);
  NumericVector FUt = FuLap(t/h, sigU);
  for (int i=0; i<nt; i++){
    R_CheckUserInterrupt();
    double ti = t[i];
    NumericVector csW = Rcpp::cos(W*ti/h);
    NumericVector snW = Rcpp::sin(W*ti/h);
    rehatFW[i] = Rcpp::sum(csW);
    imhatFW[i] = Rcpp::sum(snW);
    reYhatFW[i] = Rcpp::sum(Y*csW);
    imYhatFW[i] = Rcpp::sum(Y*snW);
  }
  for (int i=0; i<nx; i++){
    R_CheckUserInterrupt();
    NumericVector rext = Rcpp::cos(x[i]*t/h);
    NumericVector imxt = Rcpp::sin(x[i]*t/h);
    NumericVector tmp0 = rehatFW*rext + imhatFW*imxt;
    NumericVector tmp1 = -rehatFW*imxt + imhatFW*rext;
    double S0 = Rcpp::sum(tmp0*FKt/FUt)*dt/(n*h*2*PI);
    double S1 = Rcpp::sum(tmp1*FKt1/FUt)*dt/(n*h*2*PI);
    double S2 = -Rcpp::sum(tmp0*FKt2/FUt)*dt/(n*h*2*PI);
    double T0 = Rcpp::sum((reYhatFW*rext + imYhatFW*imxt)*FKt/FUt)*dt/(n*h*2*PI);
    double T1 = Rcpp::sum((-reYhatFW*imxt + imYhatFW*rext)*FKt1/FUt)*dt/(n*h*2*PI);
    res[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
  }
}

// function to estimate ghat of x using JASA when U is Gaussian
void gjasaGau(NumericVector& res, const NumericVector& x, const NumericVector& t, double dt, const NumericVector& W, 
      const NumericVector& Y, double sigU, double h){
  int nt = t.size();
  int n = W.size();
  int nx = x.size();
  NumericVector rehatFW(nt); 
  NumericVector imhatFW(nt);
  NumericVector reYhatFW(nt);
  NumericVector imYhatFW(nt);
  NumericVector FKt = FK(t);
  NumericVector FKt1= FK1(t);
  NumericVector FKt2= FK2(t);
  NumericVector FUt = FuGau(t/h, sigU);
  for (int i=0; i<nt; i++){
    R_CheckUserInterrupt();
    double ti = t[i];
    NumericVector csW = Rcpp::cos(W*ti/h);
    NumericVector snW = Rcpp::sin(W*ti/h);
    rehatFW[i] = Rcpp::sum(csW);
    imhatFW[i] = Rcpp::sum(snW);
    reYhatFW[i] = Rcpp::sum(Y*csW);
    imYhatFW[i] = Rcpp::sum(Y*snW);
  }
  for (int i=0; i<nx; i++){
    R_CheckUserInterrupt();
    NumericVector rext = Rcpp::cos(x[i]*t/h);
    NumericVector imxt = Rcpp::sin(x[i]*t/h);
    NumericVector tmp0 = rehatFW*rext + imhatFW*imxt;
    NumericVector tmp1 = -rehatFW*imxt + imhatFW*rext;
    double S0 = Rcpp::sum(tmp0*FKt/FUt)*dt/(n*h*2*PI);
    double S1 = Rcpp::sum(tmp1*FKt1/FUt)*dt/(n*h*2*PI);
    double S2 = -Rcpp::sum(tmp0*FKt2/FUt)*dt/(n*h*2*PI);
    double T0 = Rcpp::sum((reYhatFW*rext + imYhatFW*imxt)*FKt/FUt)*dt/(n*h*2*PI);
    double T1 = Rcpp::sum((-reYhatFW*imxt + imYhatFW*rext)*FKt1/FUt)*dt/(n*h*2*PI);
    res[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
  }
}

// function to estimate ghat of x using OURS when error is Laplace
void gnewLap(NumericVector& ghatofx, const NumericVector& x, const NumericVector& input, const NumericVector& output, 
        double beta, double beta2, const NumericVector& mconst, const NumericVector& Kinput, 
        const NumericVector& W, const NumericVector& Y, double sigU, double h){
  int m = input.size();
  int m_mid = m/2 +1; 
  int n = W.size();
  NumericVector gWinput(m);
  NumericVector fWinput(m);
  double nh = n*h;
  for (int i=0; i<m; ++i){
    R_CheckUserInterrupt();
    double x0 = input[i];
    NumericVector a = (x0-W)/h;
    NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    NumericVector a1 = a0*a;
    NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    gWinput[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fWinput[i] = S0;
  }
  
  // Find the support for CF of KK
  int indexl = (int)(round(-1.0/h/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h/beta2+m_mid)); 
  arma::vec support = arma::ones<vec>(m);
  for (int i=0; i<(indexl-1); ++i) {support[i]=0;}
  for (int i=indexu; i<m; ++i) {support[i]=0;}
  
  // Make arma objects
  arma::vec gWin(gWinput.begin(), m, false);
  arma::vec fWin(fWinput.begin(), m, false);
  arma::vec mcon(const_cast<NumericVector&>(mconst).begin(), m, false);
  
  // FFT for fW
  arma::vec Xin=mcon%fWin;
  arma::cx_vec FfW = beta*mcon%arma::fft( Xin )%support;
  
  // FFT for gW*fW
  arma::cx_vec FgWfW=beta*mcon%arma::fft( mcon%gWin%fWin )%support;
  
  // FfU
  NumericVector FfU=FuLap(output, sigU);
  arma::vec FfU2(FfU.begin(), FfU.size(), false);
  
  // inverse FFT to get fX
  arma::cx_vec Fratio=(FfW/FfU2%support);
  arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
  
  // inverse FFT to get gX*fX
  Fratio=(FgWfW/FfU2%support);
  arma::cx_vec gXfX=mcon/beta%arma::ifft( mcon%Fratio);
  
  // estimate of gX 
  arma::vec ghat = arma::real(gXfX)/arma::real(fXF);
  int nx = x.size();
  for (int i=0; i<nx; ++i){
    int ind = (int)(x[i]/beta+m_mid)-1;
    ghatofx[i] = ghat[ind];
  }
}

// function to estimate ghat of x using OURS when error is Gaussian
void gnewGau(NumericVector& ghatofx, const NumericVector& x, const NumericVector& input, const NumericVector& output, 
        double beta, double beta2, const NumericVector& mconst, const NumericVector& Kinput, 
        const NumericVector& W, const NumericVector& Y, double sigU, double h){
  int m = input.size();
  int m_mid = m/2 +1; 
  int n = W.size();
  NumericVector gWinput(m);
  NumericVector fWinput(m);
  double nh = n*h;
  for (int i=0; i<m; ++i){
    R_CheckUserInterrupt();
    double x0 = input[i];
    NumericVector a = (x0-W)/h;
    NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    NumericVector a1 = a0*a;
    NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    gWinput[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fWinput[i] = S0;
  }
  
  // Find the support for CF of KK
  int indexl = (int)(round(-1.0/h/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h/beta2+m_mid)); 
  
  // Make arma objects
  arma::vec gWin(gWinput.begin(), m, false);
  arma::vec fWin(fWinput.begin(), m, false);
  arma::vec mcon(const_cast<NumericVector&>(mconst).begin(), m, false);
  
  // FFT for fW
  arma::vec Xin=mcon%fWin;
  arma::cx_vec FfW = beta*mcon%arma::fft( Xin );
  
  // FFT for gW*fW
  arma::cx_vec FgWfW=beta*mcon%arma::fft( mcon%gWin%fWin );
  
  // FfU
  NumericVector FfU=FuGau(output, sigU);
  arma::vec FfU2(FfU.begin(), FfU.size(), false);
  //(FfU2.elem( arma::find_nonfinite(FfU2>1e-30) )).fill(0);
  
  // inverse FFT to get fX
  arma::cx_vec Fratio=(FfW/FfU2);
  for (int i=0; i<(indexl-1); ++i) {Fratio[i]=0;}
  for (int i=indexu; i<m; ++i) {Fratio[i]=0;}
  arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
  
  // inverse FFT to get gX*fX
  Fratio=(FgWfW/FfU2);
  for (int i=0; i<(indexl-1); ++i) {Fratio[i]=0;}
  for (int i=indexu; i<m; ++i) {Fratio[i]=0;}
  arma::cx_vec gXfX=mcon/beta%arma::ifft( mcon%Fratio);
  
  // estimate of gX 
  arma::vec ghat = arma::real(gXfX)/arma::real(fXF);
  int nx = x.size();
  for (int i=0; i<nx; ++i){
    int ind = (int)(x[i]/beta+m_mid)-1;
    ghatofx[i] = ghat[ind];
  }
}
