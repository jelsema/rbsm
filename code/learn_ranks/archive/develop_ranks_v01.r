


# Setup Environment ---------------------------------------------------


library("tidyverse")

library("rsopa")

# install.packages("Rfit")
# install.packages("npsm")

library("Rfit")
library("npsm")



rm( list=ls() )


# Attempt R-estimators  -------------------------------------------------


## One-sample -----------------

### Define functions ------------------------------------

signedrank_dispersion <- function( theta, x, sfun, pfun ){
  nn <- length(x)
  dd <- (x - theta)
  dp <- abs(dd)
  rr <- rank( dp )
  sum( sfun(rr, pfun)*dp )
}


phi_1s_wilcox <- function(u){ u }
phi_1s_sign <- function(u){ rep(1,length(u)) }
score_1s <- function( r, pfun){ pfun( r / (length(r)+1) ) }


srd_solve <- function(x, sfun=score_1s, pfun=phi_1s_wilcox ){
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


### Example of functions ------------------------------------

set.seed(42)
xx <- sort(rlnorm( 10, 1, 1 ))

data.frame(
  Intent = c("HL Est", "Median"),
  Direct = c( signedrank(xx), median(xx) ),
  Manual = c( srd_solve(xx), 
              srd_solve(xx, pfun=phi_1s_sign))
)

### Illustrations ------------------------------------











## Regression -----------------

#' Page 168-169


jaeckel_dispersion <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  sum( dd*ss )
}

jaeckel_dispersion_deriv <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  sum( x * -1*ss )
}


phi_2s_wilcox <- function(u){ sqrt(12)*(u-0.5) }
phi_2s_sign <- function(u){ sign(2*u -1) }
score_2s <- function( r, pfun){ 
  a0 <- pfun( r / (length(r)+1) )
  a1 <- a0 - mean(a0)
  a1sd <- sum(a1^2)/( length(a1)+1)
  a1 / sqrt(a1sd)
}

jdf_reg <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
  theta <- optim( par=c("theta"=lm(y~x)$coef[2] ), fn=jaeckel_dispersion, 
                  method="Brent", lower=min(x), upper=max(x),
                  sfun=sfun, pfun=pfun, y=y, x=x )$par
  return( theta )
}


# Still need to figure out tau


# Example of functions ------------------------------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

plot( yy ~ xx )
Rfit::rfit( yy ~ xx, symmetric = FALSE )
Rfit::rfit( yy ~ xx, symmetric = TRUE )



betahat <- jdf_reg( y=yy, x=xx, sfun=score_2s, pfun=phi_2s_wilcox )

dd <- yy - xx*betahat

# jdf_solve( x=dd )  # Need to rerun the 1-sample/location functions

signedrank(dd)
median(dd)

srd_solve( dd )




# Try from Hettmansperger & McKean algorithm ------------------------------------
#' Section 3.7
#' 
#' Also McKean, Terpstra, & Kloke (2009)
#' 


# Q %*% t(Q)  ## H = Q1 * Q1'


## Start: rfit -----
betahat0 <- yhat0 <- NULL


Y <- as.matrix( yy, ncol=1 )
X0 <- model.matrix( yy ~ xx )
X1 <- X0[,-1,drop=FALSE]

qrx <- qr(X0)
Q0 <- qr.Q(qrx)
Q1 <- Q0[ , -1 , drop=FALSE ]


# H <- X %*% solve( t(X)%*%X ) %*% t(X)
# H <- Q %*% t(Q)

# Initialize fit
if( is.null( betahat0 ) ){
  if( is.null(yhat0) ){
    yhat0 <- Y
  }
  betahat0 <- drop(crossprod( Q1 , yhat0))
}


disp_null <- jaeckel_dispersion( 0*betahat0, Q1, yy )
disp0     <- jaeckel_dispersion( betahat0, Q1, yy )
if( disp_null < disp0 ){ betahat0 <- 0*betahat0 }

ord <- order(yy - Q1 %*% betahat0)

## Start: jaeckle optimization -----
# arguments: Q1, Y, betahat0, sfun, pfun


beta0 <- betahat0
ysd   <- sd(y2)
ys <- y2/ysd

fit0a <- optim(beta0/ysd, jaeckel_dispersion, method = "BFGS", 
      x = xx, y = ys,  gr = jaeckel_dispersion_deriv
      )

fit1a <- optim(fit0$par * ysd, jaeckel_dispersion, method = "BFGS", 
              x = xx, y = y2,  gr = jaeckel_dispersion_deriv
)

fit1a$par
# 12.29928

## Back to rfit -----
# Do convergence testing?

betahat1 <- fit1a$par

yhat <- Q1 %*% betahat1
ehat <- Y - yhat

alphahat <- srd_solve( ehat )
ehat <- ehat - alphahat
yhat <- yhat + alphahat
bhat <- lsfit(xx, yhat, intercept = FALSE)$coef
bhat

alphahat0 <- srd_solve(yy)
bhat0 <- c(alphahat0, rep(0, length(bhat) - 1))




tau0 <- gettauF0( ehat, scores = wscores )


##### 

# cbind(
#   score_2s( seq(1,length(yy)), pfun=phi_2s_wilcox ),
#   getScores(wscores, seq_len(length(yy))/(length(yy) + 1))
# )

scrs <- getScores(wscores, seq_len(length(yy))/(length(yy) + 1))
sum( xx[ order(ehat) ] * scrs * -1 )

scrs2 <- score_2s( rank(ehat), pfun=phi_2s_wilcox )

sum( xx*scrs2*-1 )

j.disp( x=matrix(xx,ncol=1), y=yy, beta=betahat0, scrs )
jaeckel_dispersion( betahat0, xx, yy )

j.grad( x=matrix(xx,ncol=1), y=yy, beta=betahat0, scrs )
jaeckel_dispersion_deriv( betahat0, xx, yy )





Rfit::rfit( yy ~ xx, symmetric = TRUE )
Rfit::rfit( yy ~ xx, symmetric = FALSE )











# Step through Rfit method ------------------------------------


betahat0 <- NULL
yhat0 <- NULL

x <- model.matrix( yy~xx )
x1 <- as.matrix( x[,colnames(x)!='(Intercept)'] )
x1 <- as.matrix(cbind(rep(1,nrow(x1)),x1))
     xm <- x[,2,drop=FALSE]

y <- yy

qrx <- qr(x1)
Q   <-as.matrix( qr.Q(qrx) )
q1  <-Q[,1]
xq  <-as.matrix( Q[,2:qrx$rank] )


# Get initial guess
if( is.null(betahat0) ) {
  if( is.null(yhat0) ){ yhat0 <- y }
  betahat0 <- drop( crossprod(xq,yhat0) )
}

## Rfit codes
if(  disp(betahat0, xq, yy, wscores) > disp(rep(0, length(betahat0)), xq, yy, wscores)  ) {
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

yhat <- xq %*% betahat
ehat <- y - yhat
alphahat <- signedrank(ehat)
ehat <- ehat - alphahat
yhat <- yhat + alphahat

bhat <- lsfit(x, yhat, intercept=FALSE)$coef


##########################################################
## Rfit Jaeckel optimization

jaeckel( as.matrix(xq[ord,]), y[ord], betahat0, scores=wscores )


j.grad <- function (x, y, beta, scrs,...) {
  x <- as.matrix(x)
  e <- y - x %*% beta
  crossprod(x[order(e),],-1*scrs )
}

j.disp <- function (beta, x, y,scrs,...) {
  e <- y - x %*% beta
  drop(crossprod(e[order(e)],scrs))
}


scrs <- getScores( wscores, seq_len(length(y))/(length(y)+1) )
sd.y <- sd(y)
ystar <- y/sd.y

fit0 <- optim(betahat0/sd.y, 
              j.disp, method = "BFGS", 
              x = xm, y = ystar, scrs=scrs,
              gr=j.grad)

fit1 <- optim(fit0$par*sd.y, 
      j.disp, method = "BFGS", 
      x = xm, y = y, scrs=scrs,
      gr=j.grad)



fit0$par
fit0a$par

fit1$par
fit1a$par






##########################################################
## Casey's codes


# jaeckel_regression <- function( theta, x, y, sfun, pfun ){
#   nn <- length(y)
#   eta <- x*theta
#   
#   dd <- (y - eta)
#   sum( score_2s(dd0, phi_2s_wilcox)*dd0 )
# }



phi_2s_wilcox <- function(u){ sqrt(12)*(u-0.5) }
phi_2s_sign <- function(u){ sign(2*u -1) }
score_2s <- function( r, pfun){ 
  a0 <- pfun( r / (length(r)+1) )
  a1 <- a0 - mean(a0)
  a1sd <- sum(a1^2)/( length(a1)+1)
  a1 / sqrt(a1sd)
}


jdf_reg <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
  theta <- optim( par=c("theta"=lm(y~x)$coef[2] ), fn=jaeckel_regression, 
                  method="Brent", lower=min(x), upper=max(x),
                  sfun=sfun, pfun=pfun, y=y, x=x )$par
  
  return( theta )
}


betahat0 <- NULL
yhat0 <- NULL

x <- model.matrix(lm(yy~xx))
x1 <- as.matrix( x[,colnames(x)!='(Intercept)'] )
x1 <- as.matrix(cbind(rep(1,nrow(x1)),x1))

nn <- length(yy)
y <- yy

qrx <- qr(x1)
Q   <-as.matrix( qr.Q(qrx) )
q1  <-Q[,1]
xq  <-as.matrix( Q[,2:qrx$rank] )


# Get initial guess
if( is.null(betahat0) ) {
  if( is.null(yhat0) ){ yhat0 <- y }
  betahat0 <- drop( crossprod(xq,yhat0) )
}



sd.y <- sd(y)
ystar <- y/sd.y

fit0 <- optim(betahat0/sd.y, jaeckel_regression, 
              method = "BFGS", 
              x = xx, y = ystar,
      gr=jaeckel_regression_deriv )


optim(fit0$par*sd.y, jaeckel_regression_deriv, 
      method = "BFGS", 
      x = xx, y = yy, 
      gr=jaeckel_regression_deriv )








# xsq <- seq( 1.6, 1.8, length.out=1000 )
# jdf <- xsq*0
# for( ii in 1:length(xsq) ){
#   jdf[ii] <- j.disp( xsq[ii], xx, yy, scrs )
# }
# 
# plot( jdf ~ xsq )
# data.frame(tt=xsq, dd=jdf) |> 
#   filter( dd == min(dd) )






######


Rfit::rfit( yy ~ xx )


# # Checking equivalency of score functions
# ys <- seq_len(length(y))/(length(y)+1)
# 
# 
# cbind(
#   getScores(wscores, ys),
#   score_2s(ys, pfun=phi_2s_wilcox )
# )
# 
# cbind( wscores@phi( ys ),
#        phi_2s_wilcox( ys )
# )











# ## This gets CLOSE to Rfit, but not quite there.
# ## Why not?
# ## Maybe same issue as before: A non-unique solution?
# jaeckel_regression0 <- function( theta, x, y, sfun, pfun ){
#   nn <- length(y)
#   dd <- y - as.matrix(cbind(rep(1,nn),x)) %*% theta
#   
#   sum( sfun(dd, pfun)*dd )
# }
# 
# 
# 
# jdf_reg <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
#   theta <- optim( par=c("theta"=lm(y~x)$coef ), fn=jaeckel_regression0, 
#                   method="L-BFGS-B", lower=-Inf, upper=Inf,
#                   sfun=sfun, pfun=pfun, y=y, x=x )$par
#   
#   return( theta )
# }
# 
# 
# jdf_reg( yy, xx )







# Test QR decomposition ------------------------------------


Y <- rnorm( 5 )
X <- matrix( runif(15,1,10), ncol=3 )


qrx <- qr(X)
Q1 <- qr.Q( qrx )


H1 <- Q1 %*% t(Q1)
H2 <- X %*% solve( t(X)%*%X ) %*% t(X)

H1 %*% Y
H2 %*% Y





