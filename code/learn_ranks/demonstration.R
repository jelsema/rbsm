
#'
#' Attempting to build my own R-estimators
#' 



# Setup Environment ---------------------------------------------------

library("tidyverse")

# install.packages("Rfit")
# install.packages("npsm")

library("Rfit")
library("npsm")
library("microbenchmark")
library("Rcpp")


rm( list=ls() )


source("functions.r")






# One-sample -----------------

set.seed(42)
xx <- sort(rlnorm( 10, 1, 1 ))

data.frame(
  Intent = c("HL Est", "Median"),
  Direct = c( signedrank(xx), median(xx) ),
  Manual = c( srd_optim(xx), 
              srd_optim(xx, pfun=phi_1s_sign))
)





# Regression -----------------

## Simple linear regression ------------------------------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

# plot( yy ~ xx )

betahat_v1 <- jdf_optim( y=yy, x=xx, sfun=score_2s, pfun=phi_2s_wilcox, method="BFGS" )
d0 <- yy - xx*betahat_v1

rmod01 <- Rfit::rfit( yy ~ xx, symmetric = TRUE )
cmod01 <- jaeckle_coef( xx, yy )

t(data.frame(
  Rfit = rmod01$coefficients ,
  Manual2 = cmod01$coef,
  Manual = c( b0 = srd_optim( d0 ),
              b1 = betahat_v1 )
))



## Multiple linear regression ------------------------------------

set.seed(10403)

nn <- 20
pp <- 3

E <- rnorm( nn, 0.25 )
X <- matrix( runif(nn*pp,1,10), ncol=pp )
X0 <- cbind( 1, X )
B <- rnorm( pp, mean=2 )
B[3] <- 0.05
Y <- 1 + X %*% B + E

rmod01 <- Rfit::rfit( Y ~ X, symmetric = TRUE )
cmod01 <- jaeckle_coef( X, Y )

betahat_v1 <- jdf_optim( Y, X, method="BFGS" )
d0 <- yy - xx*betahat_v1


t(data.frame(
  Rfit = rmod01$coefficients ,
  Manual2 = cmod01$coef,
  Manual = c( srd_optim( d0 ),
              betahat_v1 )
))



# Drop Test? ------------------------------------


rmod02 <- Rfit::rfit( Y ~ X[,c(1,2)], symmetric = TRUE )
cmod02 <- jaeckle_coef( X[,c(1,2)], Y )

drop.test( rmod01, rmod02 )

summary(rmod01)

Bhat_rfit1 <- rmod01$coefficients

d1 <- disp( rmod01$coefficients, X0, Y, scores=wscores )
d2 <- disp( rmod02$coefficients, X0[,c(1,2,3)], Y, scores=wscores )

df1 <- length(rmod01$coefficients) - length(rmod02$coefficients)
df2 <- length(rmod01$y) - length(rmod01$coefficients)
test <- ((d2-d1)/df1)/(rmod01$tauhat/2)
pval <- 1 - pf(test, df1, df2)

c( test, pval )


jaeckel_dispersion( cmod01$coefficients, X0, Y )
jaeckel_dispersion( cmod02$coefficients, X0[,c(1,2,3)], Y )



# Scale / tau ------------------------------------

sourceCpp("ksm_tauc.cpp")


# Compare the tau from the regression estimates
dd <- cmod01$residuals
rnk <- length(cmod01$coefficients)-1

data.frame(
  Method = c("gettauF0", "gettau", "ksm_tau", "cmod"),
  Value = c(
    gettauF0( dd, p=rnk ),
    gettau(   dd, p=rnk ),
    ksm_tau(  dd, p=rnk, delta=0.8, hparm=2 ),
    cmod01$tau_est
  )
)


# Just create some data and compare tau


rm( list=ls() )
source("functions.r")
sourceCpp("ksm_tau2.cpp")
sourceCpp("ksm_tau3.cpp")

nn <- 20
ee <- rnorm(nn, 0, 1 )

data.frame(
  Method = c("Rfit_f", "Rfit_f", "CMJ_r", "CMJ_c"),
  Tau = c( gettau(ee),
           gettauF0(ee),
           ksm_tau(ee),
           ksm_tau_cpp3(
             ee, 
             p=0, delta=0.8, hparm=2, 
             ordered_scores = score_2s( 1:nn, pfun=phi_2s_wilcox ),
             ordered_scores_deriv = score_2s( 1:nn, 
                                              pfun=phi_2s_wilcox, 
                                              dpfun=dphi_2s_wilcox )
           )
  )
)


## Assess the computation time for the versions of tau


sourceCpp("ksm_tau2.cpp")
sourceCpp("ksm_tau3.cpp")


nn <- 1000
microbenchmark(
  setup = {
    xx <- rnorm(nn, 0, 1 )
  },
  Rfit_r = gettau(xx),
  Rfit_f = gettauF0(xx),
  CMJ_r = ksm_tau(xx),
  CMJ_c = ksm_tau_cpp3(xx, 
              p=0, delta=0.8, hparm=0, 
              ordered_scores = score_2s( 1:nn, pfun=phi_2s_wilcox ),
              ordered_scores_deriv = score_2s( 1:nn, 
                                               pfun=phi_2s_wilcox, 
                                               dpfun=dphi_2s_wilcox )
               ),
  times = 200
)


# Just the Fortran and C++ versions
nn <- 1000
microbenchmark(
  setup = {
    xx <- rnorm(nn, 0, 1 )
  },
  Rfit_f = gettauF0(xx),
  CMJ_c = ksm_tau_cpp3(xx, 
                       p=0, delta=0.8, hparm=0, 
                       ordered_scores = score_2s( 1:nn, pfun=phi_2s_wilcox ),
                       ordered_scores_deriv = score_2s( 1:nn, 
                                                        pfun=phi_2s_wilcox, 
                                                        dpfun=dphi_2s_wilcox )
  ),
  times = 200
)



xx <- rnorm(nn, 0, 1 )

ksm_tau_cpp(xx, 
            p=0, delta=0.8, hparm=0, 
            pfun=phi_2s_wilcox, dpfun=dphi_2s_wilcox,
            score_2s=score_2s )

ksm_tau( xx )


## Simulation on the three tau functions

nn <- 30 
ord_scr <- score_2s( 1:nn, pfun=phi_2s_wilcox )
ord_drv <- score_2s( 1:nn, 
                     pfun=phi_2s_wilcox, 
                     dpfun=dphi_2s_wilcox )


nSim <- 1000
mtau <- matrix( NA, nrow=nSim, ncol=4 )

for( ii in 1:nSim ){
  
  xx <- rnorm(nn, 0, 1 )
  mtau[ii,] <- c(
    gettau(xx),
    gettauF0(xx),
    ksm_tau(xx),
    ksm_tau_cpp(xx, 
                p=0, delta=0.8, hparm=0, 
                ordered_scores = ord_scr,
                ordered_scores_deriv = ord_drv
    )
    
  )
  
}

colnames(mtau) <- c(
  "Rfit_r",
  "Rfit_f",
  "CMJ_r",
  "CMJ_c"
)
dat_mtau <- as.data.frame(mtau)

ggplot( dat_mtau, aes(x=Rfit_f) ) + 
  geom_point( aes(y=Rfit_r), color="steelblue3" ) + 
  geom_point( aes(y=CMJ_r), color="indianred3" ) + 
  geom_point( aes(y=CMJ_c), color="thistle3" )


hist( mtau[,1] - mtau[,2] )
hist( mtau[,1] - mtau[,3] )
hist( mtau[,1] - mtau[,4] )

hist( mtau[,2] - mtau[,4] )
hist( mtau[,2] - mtau[,3] )

hist( mtau[,3] - mtau[,4] )


## Correlation / Robust R^2 ------------------------------------

summary(rmod01)




