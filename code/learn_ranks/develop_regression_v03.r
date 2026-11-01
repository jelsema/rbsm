
#'
#' Attempting to build my own R-estimators
#' 

# Setup Environment ---------------------------------------------------


library("tidyverse")

# install.packages("Rfit")
# install.packages("npsm")

library("Rfit")
library("npsm")



rm( list=ls() )


# source("functions.r")
source("code/learn_ranks/functions.r")




# One-sample -----------------

set.seed(42)
xx <- sort(rlnorm( 10, 1, 1 ))

data.frame(
  Statistic = c("HL Est", "Median"),
  Direct = c( signedrank(xx), median(xx) ),
  Manual = c( srd_optim(xx), 
              srd_optim(xx, pfun=phi_1s_sign))
)





# Regression -----------------



## Setup jdf_optim ------------------------------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

# plot( yy ~ xx )


betahat_v1 <- jdf_optim( y=yy, x=xx, sfun=score_2s, pfun=phi_2s_wilcox, method="BFGS" )

dd <- yy - xx*betahat_v1

signedrank(dd)
median(dd)

data.frame(
  Rfit = Rfit::rfit( yy ~ xx, symmetric = FALAS )$coef,
  Manual = c( b0 = srd_optim( dd ),
              b1 = betahat_v1 )
)







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



# Rfit::rfit( Y ~ X[,1] + X[,2] + X[,3], symmetric = TRUE )

# jaeckel_dispersion <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
#   dd <- y - as.matrix(x) %*% beta
#   ss <- sfun( rank(dd), pfun )
#   sum( dd*ss )
# }
# 
# jaeckel_dispersion_deriv <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
#   dd <- y - as.matrix(x) %*% beta
#   ss <- sfun( rank(dd), pfun )
#   ### sum( x * -1*ss )
#   ### drop( t(x) %*% (-1*ss) )
#   drop( crossprod( x, (-1*ss) ) )
# }
# 
# jdf_optim2 <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
#   theta <- optim( par=c("theta"=lm(y~x)$coef ), fn=jaeckel_dispersion, 
#                   method="BFGS", 
#                   sfun=sfun, pfun=pfun, y=y, x=x )$par
#   return( theta )
# }


jdf_optim2 <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox, method="Brent" ){
  theta <- optim( par=c("theta"=lm(y~x)$coef[-1] ), fn=jaeckel_dispersion, 
                  method=method,
                  sfun=sfun, pfun=pfun, y=y, x=x )$par
  return( theta )
}


betahat_v1 <- jdf_optim2( Y, X, method="BFGS" )

dd <- Y - X %*% betahat_v1
signedrank(dd)
median(dd)

data.frame(
  Rfit = Rfit::rfit( Y ~ X, symmetric = TRUE )$coef,
  Manual = c( b0 = srd_optim( dd ),
              b1 = betahat_v1 )
)



## Try from Hettmansperger & McKean algorithm ------------------------------------
#' Section 3.7
#' 
#' Also McKean, Terpstra, & Kloke (2009)
#' 


# Q %*% t(Q)  # H = XQ * XQ'


### Start: rfit -----
betahat0 <- yhat0 <- NULL

# ## Commenting this out so that I can test different setups
# Y <- as.matrix( yy, ncol=1 )
# X0 <- model.matrix( yy ~ xx )
X1 <- X0[,-1,drop=FALSE]

qrx <- qr(X0)
Q0 <- qr.Q(qrx)
XQ <- Q0[ , -1 , drop=FALSE ]


# H <- X %*% solve( t(X)%*%X ) %*% t(X)
# H <- Q %*% t(Q)

# Initialize fit
if( is.null( betahat0 ) ){
  if( is.null(yhat0) ){
    yhat0 <- Y
  }
  betahat0 <- drop(crossprod( XQ , yhat0))
}


disp_null <- jaeckel_dispersion( 0*betahat0, XQ, Y )
disp_init <- jaeckel_dispersion( betahat0, XQ, Y )
if( disp_null < disp_init ){ betahat0 <- 0*betahat0 }

ord <- order(Y - XQ %*% betahat0)


# plot( yy ~ XQ )
# jaeckel_dispersion( betahat0, x=XQ, y=yy )
# jaeckel_dispersion_deriv( betahat0, x=XQ, y=yy )



### Start: jaeckle optimization -----
# arguments: XQ, Y, betahat0, sfun, pfun
control <- list( maxit = 200, reltol = .Machine$double.eps^(3/4) )

beta0 <- betahat0
ysd   <- sd(Y)
ys <- Y/ysd

fit0a <- optim(beta0/ysd, jaeckel_dispersion, method = "BFGS", 
      x = XQ, y = ys,  gr = jaeckel_dispersion_deriv,
      control = control
      )

fit1a <- optim(fit0a$par * ysd, jaeckel_dispersion, method = "BFGS", 
              x = XQ, y = Y,  gr = jaeckel_dispersion_deriv,
              control = control
)

fit1a$par
jaeckle_coef_noint( XQ, Y, betahat0 )$par
# 12.29928
# -37.963239  15.482574   5.536306  # from Rfit

# Should do convergence a convergence assessment here?
# Rfit does some stability check

# ee <- Y - XQ %*% beta0/ysd
# crossprod(
#   XQ,
#   score_2s( ee, pfun=phi_2s_wilcox )
# )



### Back to rfit -----

betahat1 <- fit1a$par

yhat <- XQ %*% betahat1
ehat <- Y - yhat

alphahat <- srd_optim( ehat )
ehat <- ehat - alphahat
yhat <- yhat + alphahat
# Should be able to replace lsfit() with a matrix / QR-approach
bhat <- lsfit(X0, yhat, intercept = FALSE)$coef 
bhat


# 0.8308056 2.9891333 1.4407253 0.6066780 
# Rfit::rfit( Y ~ X , symmetric = TRUE )
# Rfit::rfit( Y ~ X , symmetric = FALSE )

alphahat0 <- srd_optim(Y)
bhat0 <- c(alphahat0, rep(0, length(bhat) - 1))

ehat1 <- Y - X0 %*% bhat
gettauF0( ehat1, p=ncol(X), scores = wscores )

# tau0 <- gettauF0( ehat, scores = wscores )

# mod01 <- Rfit::rfit( Y ~ X , symmetric = TRUE )
# mod01$taushat


jaeckle_coef( X, Y )$coef




### Test QR decomposition ------------------------------------


Y <- rnorm( 5 )
X <- matrix( runif(15,1,10), ncol=3 )


qrx <- qr(X)
XQ <- qr.Q( qrx )


H1 <- XQ %*% t(XQ)
H2 <- X %*% solve( t(X)%*%X ) %*% t(X)

H1 %*% Y
H2 %*% Y


XR <- qr.R(qrx)

XQ %*% XR
X


### Test QR decomposition ------------------------------------




nn <- 20
pp <- 3

E <- rnorm( nn )
X <- matrix( runif(nn*pp,1,10), ncol=pp )
X0 <- cbind( 1, X )
B <- rnorm( pp, mean=2 )
Y <- 1 + X %*% B + E


fit1 <- lm( Y ~ X )
fit2 <- lm( fit1$fitted.values ~ X )


fit1$coefficients - fit2$coefficients




b1 <- lm( Y ~ X - 1 )$coef

yhat <- X %*% b1
ehat <- Y - yhat

alphahat <- mean( ehat )
ehat <- ehat - alphahat
yhat <- yhat + alphahat
# Should be able to replace lsfit() with a matrix / QR-approach
bhat <- lsfit( X1, yhat, intercept = FALSE)$coef 

lm( Y ~ X )







