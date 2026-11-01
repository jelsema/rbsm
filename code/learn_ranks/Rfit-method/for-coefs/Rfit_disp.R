disp <- function (beta, x, y, scores) 
{
  x <- as.matrix(x)
  e <- y - x %*% beta
  r <- rank(e, ties.method = "first")/(length(e) + 1)
  getScores(scores, r) %*% e
}
