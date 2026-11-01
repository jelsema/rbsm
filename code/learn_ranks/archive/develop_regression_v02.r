
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


# One-sample -----------------

## Define functions ------------------------------------

signedrank_dispersion <- function( theta, x, sfun=score_1s, pfun=phi_1s_wilcox ){
  nn <- length(x)
  dd <- (x - theta)
  dp <- abs(dd)
  rr <- rank( dp )
  sum( sfun(rr, pfun)*dp )
}


phi_1s_wilcox <- function(u){ u }
phi_1s_sign <- function(u){ rep(1,length(u)) }
score_1s <- function( r, pfun){ pfun( r / (length(r)+1) ) }


srd_optim <- function(x, sfun=score_1s, pfun=phi_1s_wilcox ){
  theta <- optim( par=c("theta"=median(x) ), fn=signedrank_dispersion, 
         method="Brent", lower=min(x), upper=max(x),
         sfun=sfun, pfun=pfun, x=x )$par
  
  nn <- length(x)
  if( identical( pfun, phi_1s_sign) & (nn %% 2 == 0)  ){
    xsrt <- sort(x)
    xm1 <- xsrt[ (nn/2) ]
    xm2 <- xsrt[ (nn/2) + 1]
    if( (xm1 <= theta) & (theta <= xm2) ){
      theta <- (xm1+xm2)/2
    } else{
      stop("Minimizing L1 loss but estimate is not between middle values")
    }
  }

  return( theta )
}


## Example of functions ------------------------------------

set.seed(42)
xx <- sort(rlnorm( 10, 1, 1 ))

data.frame(
  Intent = c("HL Est", "Median"),
  Direct = c( signedrank(xx), median(xx) ),
  Manual = c( srd_optim(xx), 
              srd_optim(xx, pfun=phi_1s_sign))
)

## Illustrations ------------------------------------








# Regression -----------------

#' Page 168-169

## Define functions -----

jaeckel_dispersion <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  sum( dd*ss )
}

jaeckel_dispersion_deriv <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  # sum( x * -1*ss )
  drop( t(x) %*% (-1*ss) )
}


phi_2s_wilcox <- function(u){ sqrt(12)*(u-0.5) }
phi_2s_sign <- function(u){ sign(2*u -1) }
score_2s <- function( r, pfun){ 
  a0 <- pfun( r / (length(r)+1) )
  a1 <- a0 - mean(a0)
  a1sd <- sum(a1^2)/( length(a1)+1)
  a1 / sqrt(a1sd)
}

jdf_optim <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
  theta <- optim( par=c("theta"=lm(y~x)$coef[2] ), fn=jaeckel_dispersion, 
                  method="Brent", lower=min(x), upper=max(x),
                  sfun=sfun, pfun=pfun, y=y, x=x )$par
  return( theta )
}


# Still need to figure out tau


## Setup jdf_optim ------------------------------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

# plot( yy ~ xx )


betahat_v1 <- jdf_optim( y=yy, x=xx, sfun=score_2s, pfun=phi_2s_wilcox )

dd <- yy - xx*betahat_v1

signedrank(dd)
median(dd)

c( b0 = srd_optim( dd ),
   b1 = betahat_v1 )

Rfit::rfit( yy ~ xx, symmetric = TRUE )$coef



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

nn <- 20
pp <- 3

E <- rnorm( nn, 0.25 )
X <- matrix( runif(nn*pp,1,10), ncol=pp )
X0 <- cbind( 1, X )
B <- rnorm( pp, mean=2 )
Y <- 1 + X %*% B + E





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


disp_null <- jaeckel_dispersion( 0*betahat0, XQ, yy )
disp_init <- jaeckel_dispersion( betahat0, XQ, yy )
if( disp_null < disp_init ){ betahat0 <- 0*betahat0 }

ord <- order(yy - XQ %*% betahat0)


# plot( yy ~ XQ )
# jaeckel_dispersion( betahat0, x=XQ, y=yy )
# jaeckel_dispersion_deriv( betahat0, x=XQ, y=yy )



### Start: jaeckle optimization -----
# arguments: XQ, Y, betahat0, sfun, pfun


beta0 <- betahat0
ysd   <- sd(yy)
ys <- yy/ysd

fit0a <- optim(beta0/ysd, jaeckel_dispersion, method = "BFGS", 
      x = XQ, y = ys,  gr = jaeckel_dispersion_deriv
      )

fit1a <- optim(fit0a$par * ysd, jaeckel_dispersion, method = "BFGS", 
              x = XQ, y = yy,  gr = jaeckel_dispersion_deriv
)

fit1a$par
# 12.29928
# Should do convergence a convergance assessment here

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

alphahat0 <- srd_optim(yy)
bhat0 <- c(alphahat0, rep(0, length(bhat) - 1))

ehat1 <- Y - X0 %*% bhat
gettauF0( ehat1, p=1, scores = wscores )

# tau0 <- gettauF0( ehat, scores = wscores )

# mod01 <- Rfit::rfit( yy ~ xx , symmetric = TRUE )
# mod01$taushat


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







