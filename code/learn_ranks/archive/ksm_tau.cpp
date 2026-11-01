#include <Rcpp.h>
using namespace Rcpp;

// ---- helper: hstar0 ----
double hstar0_cpp(const NumericVector& abs_diff_ord,
                  const NumericVector& wt_ord_cum_nrmlz,
                  double val) {
  int n = abs_diff_ord.size();
  double out = 0.0;
  for (int i = 0; i < n; ++i) {
    if (abs_diff_ord[i] <= val) {
      if (wt_ord_cum_nrmlz[i] > out) out = wt_ord_cum_nrmlz[i];
    }
  }
  return out;
}

// ---- helper: looptau0 ----
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
  
  int max_iter = m, num_iter = 0, iError = 0, iStage = 0;
  bool finished = false;
  //double tn0 = tn_guess;
  
  if (hval == delta) {
    finished = true;
    num_iter = 1;
  }
  
  if (!finished) iStage = (hval > delta) ? -1 : 1;
  
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
    _["tndelta"] = tn_guess,
    _["hval"] = hval,
    _["iError"] = iError,
    _["num_iter"] = num_iter,
    _["finished"] = finished
  );
}

// ---- exported: ksm_tau ----
// [[Rcpp::export]]
double ksm_tau_cpp(NumericVector ehat,
                   double p = 0.0,
                   double delta = 0.8,
                   double hparm = 2.0,
                   Function pfun = R_NilValue,
                   Function dpfun = R_NilValue,
                   Function score_2s = R_NilValue) {
  
  int n = ehat.size();
  
  IntegerVector idx = seq_len(n);
  NumericVector ordered_scores = score_2s(idx, pfun);
  NumericVector ordered_scores_deriv = score_2s(idx, pfun, dpfun);
  
  double const0 = ordered_scores[n - 1] - ordered_scores[0];
  NumericVector vphi0 = ordered_scores_deriv;
  double cn0 = sum(vphi0) / const0;
  
  int npairs = n * (n - 1) / 2;
  NumericVector abs_diff(npairs), wt(npairs);
  
  int k = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      abs_diff[k] = std::abs(ehat[i] - ehat[j]);
      wt[k] = (vphi0[i] + vphi0[j]) / const0;
      k++;
    }
  }
  
  IntegerVector o = seq_len(npairs) - 1;
  std::sort(o.begin(), o.end(),
            [&](int a, int b) { return abs_diff[a] < abs_diff[b]; });
  
  NumericVector abs_diff_ord(npairs), wt_ord(npairs);
  for (int i = 0; i < npairs; ++i) {
    abs_diff_ord[i] = abs_diff[o[i]];
    wt_ord[i] = wt[o[i]];
  }
  
  NumericVector wt_ord_cum_nrmlz = cumsum(wt_ord) / ((n - 1.0) * cn0);
  
  List tn_out = looptau0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, delta);
  double tn = as<double>(tn_out["tndelta"]) / std::sqrt((double)n);
  double pn = hstar0_cpp(abs_diff_ord, wt_ord_cum_nrmlz, tn);
  
  double gamma_est = (const0 * pn) / (2.0 * tn);
  double tau_est = std::sqrt(n / (n - p)) * (1.0 / gamma_est);
  
  Function median = Environment::stats_env()["median"];
  Function mad = Environment::stats_env()["mad"];
  
  NumericVector ehat_stnd =
    (ehat - as<double>(median(ehat))) / as<double>(mad(ehat));
  
  int cnt = 0;
  for (int i = 0; i < n; ++i)
    if (std::abs(ehat_stnd[i]) < hparm) cnt++;
    
    double hubcor = std::max(cnt / (double)n, 1e-6);
    double fincor = 1.0 + (p / n) * ((1.0 - hubcor) / hubcor);
    
    return fincor * tau_est;
}
