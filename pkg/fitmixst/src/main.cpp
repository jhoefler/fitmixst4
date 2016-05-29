#ifdef INSIDE
#include <Rcpp.h>
#include <RcppGSL.h>
#include <RInside.h>                    // for the embedded R via RInside
//#include "rcpp_hello_world.h"
#include "gsl/gsl_minmax.h"
using namespace Rcpp;
using namespace std;
int main(int argc, char *argv[]) {}
    //RInside R(argc, argv);              // create an embedded R instance
    //SEXP s = rcpp_hello_world();
   /* Language call("print",s);
    call.eval();
    return 0;
    gsl_min(5,5);

}*/
#endif
