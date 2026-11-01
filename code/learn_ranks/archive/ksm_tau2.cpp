#include <Rcpp.h>
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
// hstar0
// ======================================================
double hstar0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double val) {
  
  int n = abs_diff_ord.size();
  double out = 0.0;
  
  for (int i = 0; i < n; ++i) {
    if (abs_diff_ord[i] <= val) {
      if (wt_ord_cum_nrmlz[i] > out)
        out = wt_ord_cum_nrmlz[i];
    }
  }
  
  return out;
}


// ======================================================
// looptau0
// ======================================================
List looptau0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double delta) {
  
  int m = abs_diff_ord.size();
  
  int idx_delta = std::round(delta * m);
  if (idx_delta < 1) idx_delta = 1;
  if (idx_delta > m) idx_delta = m;
  
  int idx_keep = idx_delta;
  double tn_guess = abs_diff_ord[idx_keep - 1];
  double hval = hstar0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, tn_guess);
  
  int max_iter = m;
  int num_iter = 0;
  int iError = 0;
  int iStage = 0;
  bool finished = false;
  
  if (hval == delta) {
    finished = true;
    num_iter = 1;
  }
  
  if (!finished)
    iStage = (hval > delta) ? -1 : 1;
  
  while (!finished) {
    
    num_iter++;
    
    bool condition01 =
      ((iStage == -1 && hval > delta) ||
      (iStage ==  1 && hval < delta));
    
    if (condition01) {
      
      idx_keep += iStage;
      
      if (idx_keep < 1 || idx_keep > m) {
        finished = true;
        iError = 1;
      } else {
        tn_guess = abs_diff_ord[idx_keep - 1];
        hval = hstar0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, tn_guess);
      }
      
    } else {
      finished = true;
      iError = 1;
    }
    
    if (num_iter > max_iter) {
      finished = true;
      iError = 2;
    }
  }
  
  return List::create(
    _["tndelta"]  = tn_guess,
    _["hval"]     = hval,
    _["iError"]   = iError,
    _["num_iter"] = num_iter,
    _["finished"] = finished
  );
}


// // ======================================================
// // ksm_tau (fully C++ independent)
// // ======================================================
// // [[Rcpp::export]]
// double ksm_tau_cpp(NumericVector ehat,
//                    double p,
//                    double delta,
//                    double hparm,
//                    Function pfun,
//                    Function dpfun,
//                    Function score_2s) {
//   
//   int n = ehat.size();
//   
//   IntegerVector idx(n);
//   for (int i = 0; i < n; ++i)
//     idx[i] = i + 1;
//   
//   NumericVector ordered_scores =
//     as<NumericVector>(score_2s(idx, pfun));
//   
//   NumericVector ordered_scores_deriv =
//     as<NumericVector>(score_2s(idx, pfun, dpfun));
//   
//   double const0 = ordered_scores[n - 1] - ordered_scores[0];
//   
//   double sum_deriv = 0.0;
//   for (int i = 0; i < n; ++i)
//     sum_deriv += ordered_scores_deriv[i];
//   
//   double cn0 = sum_deriv / const0;
//   
//   int npairs = n * (n - 1) / 2;
//   NumericVector abs_diff(npairs);
//   NumericVector wt(npairs);
//   
//   int k = 0;
//   for (int i = 0; i < n - 1; ++i) {
//     for (int j = i + 1; j < n; ++j) {
//       abs_diff[k] = std::abs(ehat[i] - ehat[j]);
//       wt[k] = (ordered_scores_deriv[i] +
//         ordered_scores_deriv[j]) / const0;
//       k++;
//     }
//   }
//   
//   IntegerVector ord(npairs);
//   for (int i = 0; i < npairs; ++i)
//     ord[i] = i;
//   
//   std::sort(ord.begin(), ord.end(),
//             [&](int a, int b) {
//               return abs_diff[a] < abs_diff[b];
//             });
//   
//   NumericVector abs_diff_ord(npairs);
//   NumericVector wt_ord(npairs);
//   
//   for (int i = 0; i < npairs; ++i) {
//     abs_diff_ord[i] = abs_diff[ord[i]];
//     wt_ord[i] = wt[ord[i]];
//   }
//   
//   NumericVector wt_ord_cum_nrmlz(npairs);
//   double denom = (n - 1.0) * cn0;
//   double acc = 0.0;
//   
//   for (int i = 0; i < npairs; ++i) {
//     acc += wt_ord[i];
//     wt_ord_cum_nrmlz[i] = acc / denom;
//   }
//   
//   List tn_out = looptau0_cpp(abs_diff_ord,
//                              wt_ord_cum_nrmlz,
//                              delta);
//   
//   double tn = as<double>(tn_out["tndelta"]) /
//     std::sqrt((double)n);
//   
//   double pn = hstar0_cpp(abs_diff_ord,
//                          wt_ord_cum_nrmlz,
//                          tn);
//   
//   double gamma_est = (const0 * pn) / (2.0 * tn);
//   double tau_est =
//     std::sqrt(n / (n - p)) * (1.0 / gamma_est);
//   
//   // ---- robust correction (pure C++)
//   double med = median_cpp(ehat);
//   double madv = mad_cpp(ehat, med);
//   
//   int count = 0;
//   for (int i = 0; i < n; ++i) {
//     double z = (ehat[i] - med) / madv;
//     if (std::abs(z) < hparm)
//       count++;
//   }
//   
//   double hubcor = std::max(count / (double)n, 1e-6);
//   double fincor =
//     1.0 + (p / n) * ((1.0 - hubcor) / hubcor);
//   
//   return fincor * tau_est;
// }


// ======================================================
// ksm_tau (fully C++ independent)
// ======================================================
// [[Rcpp::export]]
double ksm_tau_cpp(NumericVector ehat,
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
  
  List tn_out = looptau0_cpp(abs_diff_ord,
                             wt_ord_cum_nrmlz,
                             delta);
  
  double tn = as<double>(tn_out["tndelta"]) /
    std::sqrt((double)n);
  
  double pn = hstar0_cpp(abs_diff_ord,
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

