

#' Illustrate the way Rfit uses to estimate in the least-squares approach

nn <- 50
pp <- 3

E <- rnorm( nn, 0.25 )
X <- matrix( runif(nn*pp,1,10), ncol=pp )
X0 <- cbind( 1, X )
B <- rnorm( pp, mean=2 )
B[3] <- 0.05
Y <- 1 + X %*% B + E



qr0 <- qr(X0)
Q <- qr.Q( qr0 )
XQ <- Q[ , -1, drop=FALSE ]


bhat2 <- drop( crossprod(XQ, Y) )
yhat <- drop( XQ %*% bhat2 )


ehat <- Y - yhat
ahat <- mean(ehat)
ehat <- ehat - ahat
yhat <- yhat + ahat



bh1a <- drop( solve( t(X0)%*%X0 ) %*% (t(X0)%*%yhat ) )
bh1b <- drop( coef(lm.fit(X0,yhat)) )
bh1c <- drop( coef(lsfit(X0,yhat, intercept = FALSE)) )




t(cbind( 
  bh1a, bh1b, bh1c,
  lm(yhat ~ X)$coef
))


lm( Y ~ X )



