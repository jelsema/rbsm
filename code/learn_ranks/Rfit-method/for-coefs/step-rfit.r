


# Setup Environment ---------------------------------------------------


library("tidyverse")

library("Rfit")
library("npsm")



rm( list=ls() )




## Setup regression data ------------------------------------



### SLR ------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

Y <- as.matrix( yy, ncol=1 )
X0 <- model.matrix( yy ~ xx )
X <- X0[,2]


### MLR ------------

set.seed(10403)

nn <- 20
pp <- 3

E <- rnorm( nn, 0.25 )
X <- matrix( runif(nn*pp,1,10), ncol=pp )
X0 <- cbind( 1, X )
B <- rnorm( pp, mean=2 )
Y <- 1 + X %*% B + E









# Step through Rfit method ------------------------------------
y <- Y
x  <- X0
x1 <- X0

betahat0 <- NULL
yhat0    <- NULL

# x <- model.matrix( yy~xx )
# x1 <- as.matrix( x[,colnames(x)!='(Intercept)'] )
# x1 <- as.matrix(cbind(rep(1,nrow(x1)),x1))
# xm <- x[,2,drop=FALSE]



qrx <- qr(x1)
Q   <- as.matrix( qr.Q(qrx) )
q1  <- Q[,1]
xq  <- as.matrix( Q[,2:qrx$rank] )


# Get initial guess
if( is.null(betahat0) ) {
  if( is.null(yhat0) ){ yhat0 <- y }
  betahat0 <- drop( crossprod(xq,yhat0) )
}

## Rfit codes
if(  disp(betahat0, xq, y, wscores) > disp(rep(0, length(betahat0)), xq, y, wscores)  ) {
  betahat0 <- rep(0, length(betahat0))
}

ord <- order( y - xq%*%betahat0 )

fit <- jaeckel( as.matrix(xq[ord,]), y[ord], betahat0, scores=wscores)
# if( fit$convergence != 0 ) {
#   ord <- order( y - xq%*%fit$par )
#   fit2 <- jaeckel(as.matrix(xq[ord,]), y[ord], jitter(fit$par), scores=wscores )
#   if( fit$convergence != 0 ) {
#     warning("rfit: Convergence status not zero in jaeckel")
#     if( fit2$value < fit$value ) fit <- fit2
#   } else {
#     fit <- fit2
#   }
#   rm(fit2)
# }

betahat <- fit$par
# -37.963239  15.482574   5.536306

yhat <- xq %*% betahat
ehat <- y - yhat
alphahat <- signedrank(ehat)
ehat <- ehat - alphahat
yhat <- yhat + alphahat

bhat <- lsfit(x, yhat, intercept=FALSE)$coef
# 0.8309173 2.9891214 1.4407286 0.6066676 








## Rfit Jaeckel optimization -----

jaeckel( as.matrix(xq[ord,]), y[ord], betahat0, scores=wscores )


j.grad <- function (x, y, beta, scrs, ...) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  crossprod(x[order(e),],-1*scrs )
}

j.disp <- function (beta, x, y, scrs, ...) {
  e <- y - x %*% beta
  drop(crossprod(e[order(e)],scrs))
}


scrs <- getScores( wscores, seq_len(length(y))/(length(y)+1) )
sd.y <- sd(y)
ystar <- y/sd.y

fit0 <- optim(betahat0/sd.y, 
              j.disp, method = "BFGS", 
              x = xq[ord,], y = ystar[ord], scrs=scrs,
              gr=j.grad)
# -3.9452512  1.6089838  0.5753589

fit1 <- optim(fit0$par*sd.y, 
      j.disp, method = "BFGS", 
      x = xq[ord,], y = y[ord], scrs=scrs,
      gr=j.grad)
# -37.963357  15.482513   5.536401  # rfit
# -37.963357  15.482513   5.536401  # casey

# j.grad( xq[ord,], y, betahat0/sd.y, scrs )
# e <- y - xq[ord,] %*% betahat0/sd.y
# crossprod( xq[order(e),],-1*scrs )

# jd1 <- j.disp( betahat0/sd.y, x=xq[ord,], y = ystar[ord,], scrs=scrs )
# jd2 <- jaeckel_dispersion( beta0/ysd, XQ, ys )
# jd1 - jd2
# 
# jg1 <- j.grad( betahat0/sd.y, x=xq[ord,], y = ystar[ord,], scrs=scrs )
# jg2 <- jaeckel_dispersion_deriv( beta0/ysd, XQ, ys )
# jg1 - jg2
# 
# #### Investigate differenec in dispersion function ------
# 
# ##### rfit version
# 
# e <- ystar - xq[ord,] %*% (betahat0/sd.y)
# drop(crossprod(e[order(e)],scrs))
# 
# ##### casey version
# 
# dd <- ys - XQ %*% (beta0/ysd)
# ss <- sfun( rank(dd), pfun )
# sum( dd*ss )





fit0$par
fit0a$par

fit1$par
fit1a$par





## Incrementally running rfit -----






rfit_check <- function (formula, data, subset, yhat0 = NULL, scores = Rfit::wscores, 
                          symmetric = FALSE, TAU = "F0", betahat0 = NULL, ...){
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  
  # return(
  #   list(
  #     mf = mf,
  #     m = m,
  #     x = x
  #   )
  # )
  
  
  x1 <- as.matrix(x[, colnames(x) != "(Intercept)"])
  x1 <- as.matrix(cbind(rep(1, nrow(x1)), x1))
  y <- model.response(mf)
  qrx <- qr(x1)
  Q <- as.matrix(qr.Q(qrx))
  q1 <- Q[, 1]
  xq <- as.matrix(Q[, 2:qrx$rank])

  # return(
  #   list(
  #     x1=x1, y=y, Q=Q, q1=q1, xq=xq
  #   )
  # )

  if (is.null(betahat0)) {
    if (is.null(yhat0))
      yhat0 <- y
    betahat0 <- drop(crossprod(xq, yhat0))
  }
  
  disp0 <- disp(rep(0, length(betahat0)), xq, y, scores)
  disp1 <- disp(betahat0, xq, y, scores)

  if (disp(betahat0, xq, y, scores) > disp(rep(0, length(betahat0)), xq, y, scores)) {
    betahat0 <- rep(0, length(betahat0))
  }
  ord <- order(y - xq %*% betahat0)
  
  # return(
  #   list(
  #     yhat0=yhat0, betahat0=betahat0,
  #     disp_init=disp1, disp_null=disp0
  #   )
  # )
  
  fit <- jaeckel(as.matrix(xq[ord, ]), y[ord], betahat0, scores = scores,  ...)
  
  # return(
  #   list(
  #     fit=fit,
  #     xarg = as.matrix(xq[ord, ]),
  #     yarg = y[ord],
  #     yorg = y,
  #     bhat0arg = betahat0,
  #     ord = ord,
  #     xq = xq
  #   )
  # )
  
  betahat <- fit$par
  yhat <- xq %*% betahat
  ehat <- y - yhat
  alphahat <- ifelse(symmetric, signedrank(ehat), median(ehat))
  
  ehat <- ehat - alphahat
  yhat <- yhat + alphahat
  bhat <- lsfit(x, yhat, intercept = FALSE)$coef
  
  
  alphahat0 <- ifelse(symmetric, signedrank(y), median(y))
  bhat0 <- c(alphahat0, rep(0, length(bhat) - 1))
  
  # return(
  #   list(
  #     data.frame( yhat=yhat, ehat=ehat),
  #     alphahat=alphahat, betahat = bhat,
  #     x=x,
  #     alphahat0, bhat0
  #   )
  # )
  
  r.gettau <- switch(TAU, F0 = gettauF0, R = gettau, N = function(...) NA)
  tauhat <- r.gettau(ehat, ncol(xq), scores, ...)
  if (symmetric) {
    taushat <- tauhat
  }
  else {
    taushat <- taustar(ehat, qrx$rank)
  }
  
  return(
    list(
      data.frame( yhat=yhat, ehat=ehat),
      alphahat=alphahat, betahat = bhat,
      x=x,
      tau = taushat
    )
  )
  
}

tmp <- rfit_check( yy ~ xx, symmetric = TRUE, scores=wscores )
tmp

plot( tmp$yorg[tmp$ord] ~ tmp$xq[tmp$ord] )
plot( tmp$yarg ~ tmp$xarg )




## Incrementally running jaeckel ------


jaeckle_check <- function (x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], scores = Rfit::wscores, 
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
  
  return(
    list(
      sdy = sd.y,
      ystar = ystar,
      fot0 = fit0,
      scrs = scrs
    )
  )
  
  # optim(fit0$par * sd.y, j.disp, method = "BFGS", x = x, y = y, 
  #       scrs = scrs, gr = j.grad, control = control, ...)
}

smp <- jaeckle_check(
  x = tmp$xarg, 
  y = tmp$yarg, 
  beta0 = tmp$bhat0arg, 
  scores = Rfit::wscores
) 
smp

j.grad <- function(x, y, beta, scrs, ...) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  crossprod(x[order(e), ], -1 * scrs)
}
j.disp <- function(beta, x, y, scrs, ...) {
  e <- y - x %*% beta
  drop(crossprod(e[order(e)], scrs))
}



j.grad <- function(x, y, beta, scrs, ...) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  crossprod(x[order(e), ], -1 * scrs)
}
j.disp <- function(beta, x, y, scrs, ...) {
  e <- y - x %*% beta
  drop(crossprod(e[order(e)], scrs))
}


j.disp(
  beta=tmp$bhat0arg, x = tmp$xarg, y=tmp$yarg, scrs=smp$scrs
)
j.grad(
  x = tmp$xarg, y=tmp$yarg, beta=tmp$bhat0arg, scrs=smp$scrs
)



# Develop estimate of tau ------------------------------------





































