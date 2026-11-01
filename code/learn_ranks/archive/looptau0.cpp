#include <Rcpp.h>
using namespace Rcpp;

// assumes hstar0_cpp is already defined in this file
double hstar0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double val);

List looptau0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double delta) {
  
  int m = abs_diff_ord.size();
  int idx_delta = std::round(delta * m);   // R-style rounding
  if (idx_delta < 1) idx_delta = 1;
  if (idx_delta > m) idx_delta = m;
  
  int idx_keep = idx_delta;
  double tn_guess = abs_diff_ord[idx_keep - 1];
  double hval = hstar0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, tn_guess);
  
  int max_iter = m;
  int num_iter = 0;
  bool finished = false;
  int iError = 0;
  int iStage = 0;
  double tn0 = tn_guess;
  
  if (hval == delta) {
    tn0 = tn_guess;
    finished = true;
    num_iter = 1;
  }
  
  if (!finished) {
    if (hval > delta) {
      iStage = -1;
    } else {
      iStage = 1;
    }
  }
  
  while (!finished) {
    num_iter++;
    
    bool condition01 =
      ((iStage == -1 && hval > delta) ||
      (iStage ==  1 && hval < delta));
    
    if (condition01) {
      idx_keep += iStage;
      
      if (idx_keep < 1 || idx_keep > m) {
        tn0 = tn_guess;
        finished = true;
        iError = 1;
      } else {
        tn_guess = abs_diff_ord[idx_keep - 1];
        hval = hstar0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, tn_guess);
      }
    } else {
      tn0 = tn_guess;
      finished = true;
      iError = 1;
    }
    if (num_iter > max_iter) {
      tn0 = tn_guess;
      finished = true;
      iError = 2;
    }
  }
  return List::create(
    _["tndelta"] = tn0,
    _["hval"] = hval,
    _["iError"] = iError,
    _["num_iter"] = num_iter,
    _["finished"] = finished
  );
}

