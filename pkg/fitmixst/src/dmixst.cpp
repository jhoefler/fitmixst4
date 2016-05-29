
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP file77c6a95fb7( SEXP x, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu) ;
}

// definition

RcppExport SEXP file77c6a95fb7( SEXP x, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu ){
BEGIN_RCPP

NumericVector xx(x);

NumericVector pro2(pro);

NumericVector mu_neu(mu);
NumericVector Sigma_neu(Sigma);
NumericVector delta_neu(delta);
NumericVector nu_neu(nu);


int l=xx.size();
NumericVector yhat (l);

int m=pro2.size();
int p=1;

for(int i=0;i<l;i++){
yhat(i)=0;
for(int j=0;j<m;j++){
yhat(i)+=pro2(j)*2*Rf_dt((xx(i)-mu_neu(j))/sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j)),nu_neu(j),0)/
sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*Rf_pt(delta_neu(j)*1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*
(xx(i)-mu_neu(j)) *sqrt((nu_neu(j)+p)/(nu_neu(j)+(xx(i)-mu_neu(j))*1/(Sigma_neu(j)+delta_neu(j)* delta_neu(j))*
(xx(i)-mu_neu(j))))/sqrt(1-delta_neu(j)*1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*delta_neu(j)),nu_neu(j)+p,1,0);
};
};


return yhat;

END_RCPP
}



