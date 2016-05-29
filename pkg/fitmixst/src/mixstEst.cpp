
// includes from the plugin
#include <RcppGSL.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
#include <RcppGSL.h>
#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
//struct my_f_params {double nu; double a; double b; int m;};

using namespace Rcpp;

NumericVector dmixstc(NumericVector xx, NumericVector pro2, NumericVector mu_neu, NumericVector 
                      Sigma_neu, NumericVector delta_neu, NumericVector nu_neu){
    
    
    int l=xx.size();
    NumericVector yhat (l);
    
    int m=pro2.size();
    int p=1;
    
    
    for(int j=0;j<m;j++){
        yhat+=pro2(j)*2*dt((xx-mu_neu(j))/sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j)),nu_neu(j))/sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*pt(delta_neu(j)*1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))* (xx-mu_neu(j)) *sqrt((nu_neu(j)+p)/(nu_neu(j)+(xx-mu_neu(j))*1/(Sigma_neu(j)+delta_neu(j)* delta_neu(j))*(xx-mu_neu(j))))/sqrt(1-delta_neu(j)*1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*delta_neu(j)),nu_neu(j)+p);
    };
    
    
    return yhat;
    
}


NumericVector trunctm1(NumericVector mu, NumericVector sigma, double nu){
    double tmp = Rf_gammafn((nu+1)/2)/(Rf_gammafn(nu/2)*sqrt(PI*nu));
    NumericVector m1 = mu+(tmp*nu*sigma/(nu-1))
    /(((1-pt((0-mu)/sigma,nu)))*pow((1+pow(mu,2)/(nu*pow(sigma,2))),((nu-1)/2)));
    return m1;
}

NumericVector trunctm2(NumericVector mu, NumericVector sigma, double nu){
    double tmp = Rf_gammafn((nu+1)/2)/(Rf_gammafn(nu/2)*sqrt(PI*nu));
    double tmp2 = Rf_gammafn((nu-1)/2)/(Rf_gammafn((nu-2)/2)*sqrt(PI*(nu-2)));
    NumericVector m2 = -pow(sigma,2)*nu - pow(mu,2) + sigma*nu*tmp/tmp2*sqrt(nu/(nu-2))*sigma*
    (1-pt((0-mu)/(sqrt(nu/(nu-2))*sigma),nu-2))/(1-pt((0-mu)/sigma,nu))
    +2*mu*trunctm1(mu,sigma,nu);
    return m2;
    
}

struct my_f_params {double k;};


double nuf(double nu, void *p){
    struct my_f_params * params = (struct my_f_params *) p;
    double k =(params -> k);
    
    double y=log(nu/2)-Rf_digamma(nu/2)-k;
    return y;
}



double roots (double k)
{
    int status;
    int iter = 0, max_iter = 1000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;// r_expected = sqrt (5.0);
    double x_lo = 0.1, x_hi = 1000.0;
    gsl_function F;
    struct my_f_params params = {k};
    
    F.function = &nuf;
    F.params = &params;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    /*
     printf ("using %s method\n", 
     gsl_root_fsolver_name (s));
     
     printf ("%5s [%9s, %9s] %9s %10s %9s\n",
     "iter", "lower", "upper", "root", 
     "err", "err(est)");
     */
    
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.0001);
        
        
        
        /*printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
         iter, x_lo, x_hi,
         r, r - r_expected, 
         x_hi - x_lo);*/
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free (s);
    
    if (status != GSL_SUCCESS)
        //printf ("not Converged\n");
    
    return r;
}







/*double f (double x, void * p) {
 struct my_f_params * params = (struct my_f_params *)p;
 double nu = (params->nu);
 double a = (params->a);
 double b = (params->b);
 double m = (params->m);
 double f = pow(x,m)*Rf_dt((x-a)/b,nu,0);
 return f;
 
 }  */

// declarations
extern "C" {
    SEXP file44b370f6487( SEXP y, SEXP g, SEXP itermax, SEXP error, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose) ;
}

// definition

SEXP file44b370f6487( SEXP y, SEXP g, SEXP itermax, SEXP error, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose ){
    BEGIN_RCPP
#include <Rcpp.h>
    
    NumericVector yy(y);
    int gg=as<int>(g);
    int iitermax=as<int>(itermax);
    double eerror=as<double>(error);
    NumericVector pro_neu(pro);
    NumericVector mu_neu(mu);
    NumericVector Sigma_neu(Sigma);
    NumericVector delta_neu(delta);
    NumericVector nu_neu(nu);
    int verb=as<int>(verbose);
    
    
    
    int p=1;
    int n=yy.size();
    double lik = 1.0;
    double lik_neu = 1.0;
    
    
    int iter = 0;
    
    NumericVector Omega(gg), Omega_inv(gg), Lambda(gg), q(n), d(n), y_star(n), 
    e1(n), e2(n), e3(n), e4(n) ,t1(n), t2(n), t3(n), tmp(n), S2(n), S3(n), pro2(n);
    
    
    
    while(iter < iitermax && fabs(lik / lik_neu -1) > eerror || iter==0)
    {
        NumericVector pro=pro_neu;
        NumericVector mu=mu_neu;
        NumericVector delta=delta_neu;
        NumericVector Sigma=Sigma_neu;
        NumericVector nu=nu_neu;
        
        NumericVector mu1(gg);
        NumericVector Sigma1(gg);
        NumericVector delta1(gg);
        NumericVector nu1(gg);
        NumericVector pro1(gg);
        
        
        for(int i=0; i<gg; i++){
            
            Omega(i)=Sigma(i)+delta(i)*delta(i);
            Omega_inv(i)=1/Omega(i);
            
            Lambda(i)=1-delta(i)*Omega_inv(i)*delta(i);
            
            pro2=dmixstc(yy,wrap(pro(i)), wrap(mu(i)), wrap(Sigma(i)), wrap(delta(i)), wrap(nu(i)))/dmixstc(yy, pro, mu, Sigma, delta, nu);
            
            q=delta(i)*Omega_inv(i)*(yy-mu(i));
            d=(yy-mu(i))*Omega_inv(i)*(yy-mu(i));
            y_star=q*sqrt((nu(i)+p)/(nu(i)+d));
            
            t1=pt(q*sqrt((nu(i)+p+2)/(nu(i)+d))/sqrt(Lambda(i)),nu(i)+p+2);
            t2=pt(y_star/sqrt(Lambda(i)),nu(i)+p);
            
            
            e2=(nu(i)+p)/(nu(i)+d)*t1/t2;
            e1=(e2 - log((nu(i)+d)/2) - (nu(i)+p)/(nu(i)+d)+Rf_digamma((nu(i)+p)/2));
            
            tmp=(nu(i)+p)/(nu(i)+d)*pt(q/sqrt((nu(i)+d)/(nu(i)+p+2)*Lambda(i)),nu(i)+p+2)/pt(y_star/sqrt(Lambda(i)),nu(i)+p);
            
            //moments of truncated t distribution
            S2=trunctm1(q,sqrt(((nu(i)+d)/(nu(i)+p+2))*Lambda(i)),nu(i)+p+2);
            S3=trunctm2(q,sqrt(((nu(i)+d)/(nu(i)+p+2))*Lambda(i)),nu(i)+p+2);
            
            e3=tmp*S2;
            e4=tmp*S3;
            
            
            //updating parameter
            
            mu1(i)=sum(pro2*(e2*yy-delta(i)*e3))/sum(pro2*e2);
            delta1(i)=sum(pro2*(yy-mu1(i))*e3)/sum(pro2*e4);
            Sigma1(i)=1/sum(pro2)*sum(pro2*(delta1(i)*e4*delta1(i)-(yy-mu1(i))*e3*delta1(i)+(yy-mu1(i))*(yy-mu1(i))*e2-delta1(i)*e3*(yy-mu1(i))));
            
            double k= 1/sum(pro2)*sum(pro2*(e2-e1-1));
            
            nu1(i)=roots(k);                          
            
            pro1(i)=sum(pro2)/n;
            
            
            
        };
        
        mu_neu=mu1;
        delta_neu=delta1;
        Sigma_neu=Sigma1;
        nu_neu=nu1;
        pro_neu=pro1;
        
        
        

        lik = lik_neu;
        
        lik_neu = sum(log(dmixstc(yy,pro_neu,mu_neu,Sigma_neu,delta_neu,nu_neu)));
        
        iter++;
    }
    Rcpp::List coefs;
    coefs=Rcpp::List::create(Rcpp::Named("pro") = pro_neu,
                              Rcpp::Named("mu") = mu_neu,
                              Rcpp::Named("Sigma") = Sigma_neu,
                              Rcpp::Named("delta") = delta_neu,
                              Rcpp::Named("nu") = nu_neu);
    
    return Rcpp::List::create(Rcpp::Named("coefficients") = coefs,
                              Rcpp::Named("logLik") = lik_neu,
                              Rcpp::Named("lik") = lik,
                              Rcpp::Named("S2") = S2,
                              Rcpp::Named("S3") = S3,
                              Rcpp::Named("iter") = iter);
    
    END_RCPP
}



