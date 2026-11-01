#include <Rcpp.h>
//#include <algorithm>
using namespace Rcpp;

// ======================================================
// Utility: median (no full sort)
// ======================================================
double median_cpp(NumericVector x) {
  
  int n = x.size();
  NumericVector tmp = clone(x);
  
  int mid = n / 2;
  
  std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
  double med = tmp[mid];
  
  if (n % 2 == 0) {
    std::nth_element(tmp.begin(), tmp.begin() + mid - 1, tmp.end());
    med = 0.5 * (med + tmp[mid - 1]);
  }
  
  return med;
}

// ======================================================
// Utility: MAD (consistent with default R mad(), constant=1)
// ======================================================
double mad_cpp(NumericVector x, double med) {
  
  int n = x.size();
  NumericVector dev(n);
  
  for (int i = 0; i < n; ++i)
    dev[i] = std::abs(x[i] - med);
  
  double mad = median_cpp(dev);
  
  return mad;
}



// ======================================================
// hstar
// ======================================================
// [[Rcpp::export]]
double hstar_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double val) {
  
  int m = abs_diff_ord.size();
  // ---- early exits
  if (val < abs_diff_ord[0])
    return 0.0;
  
  if (val >= abs_diff_ord[m - 1])
    return wt_ord_cum_nrmlz[m - 1];  // typically 1.0
  
  auto it = std::upper_bound(abs_diff_ord.begin(),
                             abs_diff_ord.end(),
                             val);
  
  int idx = (it - abs_diff_ord.begin()) - 1;
  return wt_ord_cum_nrmlz[idx];
}


// ======================================================
// looptau
// ======================================================
double looptau_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double delta) {
  
  const int m = wt_ord_cum_nrmlz.size();
  
  // ---- edge cases
  if (delta <= wt_ord_cum_nrmlz[0])
    return abs_diff_ord[0];
  
  if (delta >= wt_ord_cum_nrmlz[m - 1])
    return abs_diff_ord[m - 1];
  
  // ---- binary search for first F(x) >= delta
  auto it = std::lower_bound(wt_ord_cum_nrmlz.begin(),
                             wt_ord_cum_nrmlz.end(),
                             delta);
  
  int idx = it - wt_ord_cum_nrmlz.begin();
  // int idx = static_cast<int>(it - wt_ord_cum_nrmlz.begin());
  // return abs_diff_ord[idx];  
  double tn_val = abs_diff_ord[idx];
  return tn_val;  
  
}



// ======================================================
// ksm_tau (fully C++ independent)
// ======================================================
// [[Rcpp::export]]
double ksm_tau_cpp3(NumericVector ehat,
                   double p,
                   double delta,
                   double hparm,
                   NumericVector ordered_scores,
                   NumericVector ordered_scores_deriv
                   ) {
  
  int n = ehat.size();
  
  double const0 = ordered_scores[n - 1] - ordered_scores[0];
  
  double sum_deriv = 0.0;
  for (int i = 0; i < n; ++i)
    sum_deriv += ordered_scores_deriv[i];
  
  double cn0 = sum_deriv / const0;
  
  int npairs = n * (n - 1) / 2;
  NumericVector abs_diff(npairs);
  NumericVector wt(npairs);
  
  int k = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      abs_diff[k] = std::abs(ehat[i] - ehat[j]);
      wt[k] = (ordered_scores_deriv[i] +
        ordered_scores_deriv[j]) / const0;
      k++;
    }
  }
  
  IntegerVector ord(npairs);
  for (int i = 0; i < npairs; ++i)
    ord[i] = i;
  
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) {
              return abs_diff[a] < abs_diff[b];
            });
  
  NumericVector abs_diff_ord(npairs);
  NumericVector wt_ord(npairs);
  
  for (int i = 0; i < npairs; ++i) {
    abs_diff_ord[i] = abs_diff[ord[i]];
    wt_ord[i] = wt[ord[i]];
  }
  
  NumericVector wt_ord_cum_nrmlz(npairs);
  double denom = (n - 1.0) * cn0;
  double acc = 0.0;
  
  for (int i = 0; i < npairs; ++i) {
    acc += wt_ord[i];
    wt_ord_cum_nrmlz[i] = acc / denom;
  }
  
  double tn_raw = looptau_cpp(abs_diff_ord, wt_ord_cum_nrmlz, delta);
  double tn = tn_raw / std::sqrt( (double)n );
  
  double pn = hstar_cpp(abs_diff_ord,
                         wt_ord_cum_nrmlz,
                         tn);
  
  double gamma_est = (const0 * pn) / (2.0 * tn);
  double tau_est =
    std::sqrt(n / (n - p)) * (1.0 / gamma_est);
  
  // ---- robust correction (pure C++)
  double med = median_cpp(ehat);
  double madv = mad_cpp(ehat, med);
  
  int count = 0;
  for (int i = 0; i < n; ++i) {
    double z = (ehat[i] - med) / madv;
    if (std::abs(z) < hparm)
      count++;
  }
  
  double hubcor = std::max(count / (double)n, 1e-6);
  double fincor =
    1.0 + (p / n) * ((1.0 - hubcor) / hubcor);
  
  return fincor * tau_est;
}
















