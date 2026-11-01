
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




# Attempting refactor to speed up the C++ ---------------------

## H_of_t ------------------------------------

sourceCpp("ksm_tau3.cpp")

nn <- 20
# p=0, delta=0.8, hparm=0, 

v1 <- score_2s( 1:nn, pfun=phi_2s_wilcox )
v2 <- score_2s( 1:nn,
                pfun=phi_2s_wilcox,
                dpfun=dphi_2s_wilcox )



ee <- rnorm(nn, 0, 1 )

Rfit::gettauF0( ee, p=0, delta=0.8, hparm=2 )

ksm_tau_cpp3(ee, 
             p=0, delta=0.8, hparm=2, 
             ordered_scores = score_2s( 1:nn, pfun=phi_2s_wilcox ),
             ordered_scores_deriv = score_2s( 1:nn, 
                                              pfun=phi_2s_wilcox, 
                                              dpfun=dphi_2s_wilcox )
)

ksm_tau_fast(ee, 
             p=0, delta=0.8, hparm=2, 
             score_2s( 1:nn, pfun=phi_2s_wilcox ),
             score_2s( 1:nn, 
                       pfun=phi_2s_wilcox,
                       dpfun=dphi_2s_wilcox )
)



## Check timing ------------------------------------

# Compare the tau from the regression estimates

nn <- 200
microbenchmark(
  setup = {
    ee <- rnorm(nn, 0, 1 )
  },
  Rfit_r = Rfit::gettau( ee, p=0, delta=0.8, hparm=2 ),
  Rfit_f = Rfit::gettauF0( ee, p=0, delta=0.8, hparm=2 ),
  CMJ_r = ksm_tau_cpp3(ee, 
                       p=0, delta=0.8, hparm=2, 
                       ordered_scores = score_2s( 1:nn, pfun=phi_2s_wilcox ),
                       ordered_scores_deriv = score_2s( 1:nn, 
                                                        pfun=phi_2s_wilcox, 
                                                        dpfun=dphi_2s_wilcox )
  )
  ,
  CMJ_c = ksm_tau_fast(ee, 
                       p=0, delta=0.8, hparm=2, 
                       score_2s( 1:nn, pfun=phi_2s_wilcox ),
                       score_2s( 1:nn, 
                                 pfun=phi_2s_wilcox,
                                 dpfun=dphi_2s_wilcox )
  )
  ,
  times = 200
)
















