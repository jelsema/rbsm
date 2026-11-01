#include <Rcpp.h>
using namespace Rcpp;

double hstar0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double val) {
  int n = abs_diff_ord.size();
  double hstar_val = 0.0;
  
  for (int i = 0; i < n; ++i) {
    if (abs_diff_ord[i] <= val) {
      if (wt_ord_cum_nrmlz[i] > hstar_val) {
        hstar_val = wt_ord_cum_nrmlz[i];
      }
    }
  }
  
  return hstar_val;
}
