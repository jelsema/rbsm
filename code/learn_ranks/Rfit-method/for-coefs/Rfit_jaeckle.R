jaeckle <- function (x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], scores = Rfit::wscores, 
          control = NULL, ...) 
{
  if (is.null(control)) {
    control <- list(maxit = 200, reltol = .Machine$double.eps^(3/4))
  }
  else {
    if (is.null(control$reltol)) {
      control$reltol <- .Machine$double.eps^(3/4)
    }
    if (is.null(control$maxit)) {
      control$maxit <- 200
    }
  }
  scrs <- getScores(scores, seq_len(length(y))/(length(y) + 1))
  j.grad <- function(x, y, beta, scrs, ...) {
    x <- as.matrix(x)
    e <- y - x %*% beta
    crossprod(x[order(e), ], -1 * scrs)
  }
  j.disp <- function(beta, x, y, scrs, ...) {
    e <- y - x %*% beta
    drop(crossprod(e[order(e)], scrs))
  }
  sd.y <- sd(y)
  ystar <- y/sd.y
  fit0 <- optim(beta0/sd.y, j.disp, method = "BFGS", x = x, 
                y = ystar, scrs = scrs, gr = j.grad, control = control, 
                ...)
  optim(fit0$par * sd.y, j.disp, method = "BFGS", x = x, y = y, 
        scrs = scrs, gr = j.grad, control = control, ...)
}



