gettau <- function (ehat, p = 0, scores = Rfit::wscores, delta = 0.8, hparm = 2, 
          ...) 
{
  n = length(ehat)
  asc = getScores(scores, 1:n/(n + 1))
  ascpr = getScoresDeriv(scores, 1:n/(n + 1))
  temph = hstarreadyscr(ehat, asc, ascpr)
  abdord = temph$absdifford
  wtord = temph$wtsord
  const = temph$cn
  templ = looptau(delta, abdord, wtord, const, n)
  tn = templ$quan/sqrt(n)
  pn = hstar(abdord, wtord, const, n, tn)
  tauscr = ((asc[n] - asc[1]) * pn)/(2 * tn)
  tauscr = sqrt(n/(n - p)) * (1/tauscr)
  w = rep(0, n)
  stan = (ehat - median(ehat))/mad(ehat)
  w[abs(stan) < hparm] = 1
  hubcor = sum(w)/n
  if (hubcor < 1e-06) {
    hubcor = 1e-06
  }
  fincor = 1 + (((p)/n) * ((1 - hubcor)/hubcor))
  tauscr = fincor * tauscr
  tauscr
}










