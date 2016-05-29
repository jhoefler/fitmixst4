#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector TPCc(NumericVector thres, NumericVector y, NumericVector yhat) {
  
  int n = thres.size();  
  IntegerVector res(n);
  
  
  //NumericVector thresholds = clone(yhat);
  //std::sort(thresholds.begin(), thresholds.end());
  for( int i=0; i<n; i++) {
      res[i] = sum( (yhat > thres[i]) & (y == 1));
  }
    return res;
}  

// [[Rcpp::export]]
IntegerVector FPCc(NumericVector thres, NumericVector y, NumericVector yhat) {
  
  int n = thres.size();  
  IntegerVector res(n);
  
  
  //NumericVector thresholds = clone(yhat);
  //std::sort(thresholds.begin(), thresholds.end());
  for( int i=0; i<n; i++) {
      res[i] = sum( (yhat > thres[i]) & (y == 0));
  }
    return res;
} 
 
// [[Rcpp::export]]
int Nrelc2(NumericVector y) {
  
  int n = y.size();  
  int res = 0;
  

  for( int i=0; i<n; i++) 
  {
    for( int j=0; j<n; j++) 
    {
      res += y[i] != y[j];
    }
  }
    return res;
}

// [[Rcpp::export]]
int Nrelc1(NumericVector y) {
  
  int res;
  
  res = 2 * sum(y==1) * sum(y==0);

    return res;
}

// [[Rcpp::export]]
int Nconc(NumericVector y, NumericVector yhat) {
  
  int n = y.size();  
  int res = 0;
  

  for( int i=0; i<n; i++) 
  {
    for( int j=0; j<n; j++) 
    {
      res += (y[i] < y[j]) * (yhat[i] < yhat[j]) +
      (y[i] > y[j]) * (yhat[i] >= yhat[j]);
    }
  }
    return res;
}

// [[Rcpp::export]]
double aucc(NumericVector y, NumericVector yhat) {
  
 return 1 - (double)Nconc(y=y, yhat=yhat)  / Nrelc1(y=y);
}

// [[Rcpp::export]]
NumericMatrix netbenc(NumericVector y, NumericVector yhat) 
{
// sort and drop duplicates
NumericVector thres = yhat;
std::sort(thres.begin(), thres.end());
thres.erase( std::unique( thres.begin(), thres.end() ), thres.end() );
int n = yhat.size(); 
int ngrid = thres.size();
NumericMatrix res(ngrid, 3);
res( _,0) = thres;

NumericVector tpc = (NumericVector)TPCc(thres, y, yhat);
NumericVector fpc = n - tpc;
NumericVector cterm = thres / (1 - thres);
res( _,1)  = tpc / n  - fpc / n * cterm ;

// treating everyone
double prev = (double)sum(y==1)/n;
res( _,2)  = prev - (1-prev) * cterm ;
return res;
}
