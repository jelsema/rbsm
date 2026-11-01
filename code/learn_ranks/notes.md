

# Tasks


## Figure out methods

- Drop-in-Dispersion test
- Figure out the cluster-correlated / random effects model
- Survival analysis




## General 

- More intelligent implementation of the C++ version of the tau calculations

- Come up with a better way to handle the score functions
  - Maybe a combination of a selector and a validator
  - Selector similar to getScores
    - Put in name
    - Get out two functions, pfun and dpfun
  - Validator to allow custom scores?
    - Or not needed?


## Improvements

- Make more full-featured linear model function
  - [] Make test of global null
    - [] Get the dispersion for full and null models
  - [] Compute correlation
  - [] Compute partial tests




# Prompts to chatGPT



I'm writing an R function. My goal is to re-implement something that is already 
implemented, for the purpose of having my own implementation, and to understand 
the construction of it. 
I'm first going to describe the goal at a high level.
Then I will provide some more description of the algorithm.
Then I will provide the three C++ functions which I currently have. These functions 
get me the correct value. What I am interested in is making the computation as 
fast as possible.




The goal is to compute a robust estimate of uncertainty, like an alternative 
to the standard deviation. The current implementation is a set of 3 functions.
I do not care if these functions remain separate, or if they get integrated 
into a single function. My goal is to make it as computationally efficient as 
possible. To accomplish this, I'm writing them in C++, and using Rcpp to call 
them from R.


The "main" function is `ksm_tau`. This takes a vector of deviations called `ehat`, 
as well as two vectors ordered_scores and ordered_scores_deriv. These values 
correspond to the entries of `ehat` if that vector was sorted. There are also two 
tuning parameters, `delta` and `hparm`. The length of `ehat` is `n`.

The current implementation takes the absolute value of the pairwise differences 
in `ehat`. Using R notation, that would be `abs( ehat[i] - ehat[j])`. The task 
ultimately requires evaluating the empirical CDF (ECDF) of these pairwise deviations.
The function to evaluate the ECDF is called `hstar0`. We need the value `tn` such that 
`hstar0(tn)` evalutes to 0.8. Then we evaluate `hstar0` at the value `tn/sqrt(n)`.

There are a few additional calculations from there, but those things are the main 
computational bottlenecks, as far as I can see.






The goal of this is to come up with a robust estimate of scale. The overall process is:

- The "main" function is called `ksm_tau`
- This function takes in a vector of deviations, called `ehat`, as well as two vectors `ordered_scores` and `ordered_scores_deriv`. These values correspond to the entries of `ehat` if that vector was sorted. There are also two tuning parameters, `delta` and `hparm`. The length of `ehat` is `n`.
- The vector `ordered_scores_deriv` is what controls the weight assigned to each value of `ehat`, but it does so indirectly.
  - Obtain the maximum difference in `ordered_scores` and call this `const0`.
  - Obtain the sum of `ordered_scores_deriv`, divide by `const0`, and call the result `cn0`.
  - When obtaining the absolute value of the pairwise difference of `ehat`, the calculations are:
    `abs_diff[k] = std::abs( ehat[i] - ehat[j] )` and  `wt[k] = (ordered_scores_deriv[i] + ordered_scores_deriv[j])/const0`,
    where `k` starts at 0 and iterates as either `i` or `j` change. The vectors are defined with appropriate size in advance.
  - Then both `abs_diff` and `wt` are ordered according to `abs_diff` to obtain `abs_diff_ord` and `wt_ord`.
  - Finally, we get `wt_ord_cum_normlz` by taking the cumulative sum of `wt_ord` and dividing it by `(n-1)*cn0`.
  - The vector `wt_ord_cum_normlz` is the weight for each absolute pairwise difference.

- My current implementation.
  - I obtain the absolute value for all the pairwise differences, `abs_diff`, as well as the weights `wt`.
  - I then sort these vectors according to the magnitude of `abs_diff`, to obtain `abs_diff_ord` and `wt_ord`. 
    - To keep it separate and distinct, length of `abs_diff_ord` is `m`.
  - A function `hstar0` takes `abs_diff_ord` and the corresponding ordered weights after normalization, `wt_ord_cum_normlz` and evaluates the ECDF as a specified value given by a third argument `val`.
  - A function `looptau0` takes `abs_diff_ord` and `wt_ord_cum_normlz`, as well as the argument `delta`, to find the `delta` quantile of the ECDF. Normally `delta` will be 0.8 or 0.9. This `looptau0` calls `hstar0` several times. The output from `looptau0` is the value `tn`.
  - Then `hstar0` is evaluated at the value `tn/sqrt(n)`, where `n` is the length of `ehat`.
  - There are a few other computations, but I do not think those are significant bottlenecks.

Taking all that into consideration:

- Is it still advisable to not compute `abs_diff_ord`?
- If so, write the function or functions.








I misspoke earlier. The vectors `ordered_scores` and `ordered_scores_deriv` are not directly associated to the input `ehat`. They would be associated with that only if `ehat` was sorted. For computational reasons, these vectors are computed in advance assuming an ordered vector, since the existing implementation will later order the absolute value of the differences of `ehat`. If needed, they can be redone to be directly associated with the entries in `ehat`.

Then regarding the weights, they are not a simple multiplication. The values of `ordered_scores_deriv` start off as a set of constants (often all the same, but not always) which are not normalized. The weights get processed as:

- Obtain the maximum difference in `ordered_scores` and call this `const0`.
- Obtain the sum of `ordered_scores_deriv`, divide by `const0`, and call the result `cn0`.
- When obtaining the absolute value of the pairwise difference of `ehat`, the calculations are:
  `abs_diff[k] = std::abs( ehat[i] - ehat[j] )` and  `wt[k] = (ordered_scores_deriv[i] + ordered_scores_deriv[j])/const0`,
  where `k` starts at 0 and iterates as either `i` or `j` change. The vectors are defined with appropriate size in advance.
- Then both `abs_diff` and `wt` are ordered according to `abs_diff` to obtain `abs_diff_ord` and `wt_ord`.
- Finally, we get `wt_ord_cum_normlz` by taking the cumulative sum of `wt_ord` and dividing it by `(n-1)*cn0`.





double hstar0(const NumericVector& abs_diff_ord,
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


double looptau0(const NumericVector& abs_diff_ord,
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




double ksm_tau(NumericVector ehat,
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
    
  double tn_raw = looptau0(abs_diff_ord, wt_ord_cum_nrmlz, delta);
  double tn = tn_raw / std::sqrt( (double)n );
  
  double pn = hstar0(abs_diff_ord,
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
