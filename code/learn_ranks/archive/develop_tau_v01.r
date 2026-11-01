

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


source("functions.r")



# Regression -----------------


## Setup regression data ------------------------------------

set.seed(10403)

ee <- rlnorm( 10, 1, 1 )
ee <- ee - median(ee)
xx <- seq(1,10,1)
yy <- 3 + 1.5*xx + ee

plot( yy ~ xx )



betahat_v1 <- jdf_optim( y=yy, x=xx, sfun=score_2s, pfun=phi_2s_wilcox, method="BFGS" )

dd <- yy - xx*betahat_v1

signedrank(dd)
median(dd)

rmod01 <- Rfit::rfit( yy ~ xx, symmetric = TRUE )
cmod01 <- jaeckle_coef( xx, yy )

t(data.frame(
  Rfit = rmod01$coefficients ,
  Manual2 = cmod01$coef,
  Manual = c( b0 = srd_optim( dd ),
              b1 = betahat_v1 )
))

rr <- cmod01$residuals
cbind( rmod01$residuals, rr, rmod01$residuals - rr )

rmod01$taushat
gettauF0( rmod01$residuals, p=1 )


data.frame(
  Method = c("gettauF0", "gettau", "ksm_tau", "cmod"),
  Value = c(
    gettauF0( rr, p=1 ),
    gettau(   rr, p=1 ),
    ksm_tau(  rr, p=1, delta=0.8, hparm=2 ),
    cmod01$tau_est
  )
)

xtest <- rnorm( 20, 0 ,1 )
data.frame(
  Method = c("gettauF0", "gettau", "ksm_tau"),
  Value = c(
    gettauF0( xtest, p=0 ),
    gettau(   xtest, p=0 ),
    ksm_tau(  xtest, p=0, delta=0.8, hparm=2 )
  )
)






# Develop tau ------------------------------------



## START: gettau() -------------------
## args: ehat, p = 0, scores = Rfit::wscores, delta = 0.8, hparm = 2

# Set arguments for gettau
ehat <- ee
p <- 0
scores <- wscores
delta <- 0.8
hparm <- 2




n <- length(ee)

asc    <- getScores(scores, 1:n/(n + 1))         # just gets generic scores
ascpr  <- getScoresDeriv(scores, 1:n/(n + 1))    # just gets generic score derivatives

# START: hstarreadyscr() --------------------------------
# temph  <- hstarreadyscr(ehat, asc, ascpr)
#   abdord <- temph$absdifford
#   wtord  <- temph$wtsord
#   const  <- temph$cn

# phi_2s_wilcox(ehat)
# sort( score_2s(1:n, phi_2s_wilcox ) )
ordered_scores <- score_2s( 1:n, phi_2s_wilcox )
ordered_scores_deriv <- score_2s( 1:n, phi_2s_wilcox, dphi_2s_wilcox )
const0 <- ordered_scores[n] - ordered_scores[1]

vphi0 <- ordered_scores_deriv
cn0 <- sum(ordered_scores_deriv)/const0

pair_idx <- combn( n,2 )
vad0 <- cbind( ehat[pair_idx[1,]], ehat[pair_idx[2,]])
vwt0 <- cbind( vphi0[pair_idx[1,]], vphi0[pair_idx[2,]])/const0
vad02 <- abs(vad0[,1] - vad0[,2])
vwt02 <- vwt0[,1] + vwt0[,2]
iord <- order(vad02)
abs_diff_ord <- vad02[iord]
wt_ord <- vwt02[iord]


# I believe that these calculations are from Eqn 3.7.7 and/or 3.7.5
#  abs_diff_ord == |e(i) - e(j)|
#  wt_ord is    == ..... not sure yet


  # orde = sort(ehat)
  # const = asc[n] - asc[1]
  # # n = length(ehat)
  # # ic = (1:n)/(n + 1)
  # vphi = rep(0, n)
  # cn = 0
  # for (i in 1:n) {
  #   vphi[i] = ascpr[i]
  #   cn = cn + ascpr[i]/const
  # }
  # vad = pairup(orde)
  # vwt = pairup(vphi)/const
  # vad2 = abs(vad[, 1] - vad[, 2])
  # vwt2 = vwt[, 1] + vwt[, 2]
  # iord = order(vad2)
  # absdifford = vad2[iord]
  # wtsord = vwt2[iord]
  # list(absdifford = absdifford, wtsord = wtsord, cn = cn)


# Check with below
# cbind( abs_diff_ord, absdifford, abs_diff_ord - absdifford )
# cbind( wt_ord, wtsord, wt_ord - wtsord )




# START: looptau() --------------------------------
# Trying to estimate t_{n,delta}
templ  <- looptau(delta, abdord, wtord, const, n)
tn     <- templ$quan/sqrt(n)


m <- length(abdord) # Number of pairs ... could choose(n,2) work? Maybe some small edge case.
const <- cn
wt_ord_cum_nrmlz <- cumsum(wt_ord) / ((n - 1) * cn)

# Demonstrate my hstar() function gets the same result
# Should make a more general demonstration
# hstar0( abs_diff_ord, cumsum(wt_ord) / ((n - 1) * cn), cn, n, 0.75 )
# hstar(absdifford, wtsord, cn, n, 0.75 )

idx_delta <- round(delta*m) # replace: inddelta
tau_guess <- abdord[idx_delta]  # replace: guess
hval <- hstar0( abs_diff_ord, wt_ord_cum_nrmlz, tau_guess ) # replace: y

max_iter <- m         # replace: itmax
num_iter <- 0         # replace: ?icmax
idx_keep <- idx_delta # replace: ikeep
finished <- FALSE     # replace: gohome



if( hval == delta ){
  tau0 <- tau_guess
  finished <- TRUE
  iterations <- 1
}
if( finished == FALSE ){
  if( hval > delta ){
    iStage <- -1  # move to the left in ordered set
  } else{
    iStage <- 1   # move to the right in ordered set
  }
}
while( finished == FALSE ){
  # WoRKING ON THIS PART
  # (iStage == -1) & (hval > delta) ---> move left
  # (iStage == -1) & (hval < delta) ---> move left
  num_iter <- num_iter + 1
  condition01 <- (((iStage == -1) & (hval > delta)) |
                  ((iStage ==  1) & (hval < delta))   )
  
  
  if( condition01 ){
    idx_keep <- idx_keep + iStage  # move left or right according to stage
    if( (idx_keep<1) | (idx_keep>m) ){
      # End loop if outside of the bounds
      tau0      <- tau_guess
      finished <- TRUE
      iError   <- 1
    } else{
      tau_guess = abdord[idx_keep]
      hval = hstar0(abs_diff_ord, wt_ord_cum_nrmlz, tau_guess)
    }
  } else{
    # End loop if:
    # - moving left but hval < delta
    # - moving right but hval > delta
    tau0      <- tau_guess
    finished <- TRUE
    iError   <- 1
  }
  
  if( num_iter > max_iter ){
    tau0      <- tau_guess
    finished <- TRUE
    iError   <- 2
  }
}


# Demonstrate my looptau() function gets the same result
# Should make a more general demonstration
# looptau0( abs_diff_ord, cumsum(wt_ord) / ((n - 1) * cn), delta=0.79 )
# looptau(0.79, abdord, wtord, const, n)


xtest <- rnorm( 8, 0, 1 )

ksm_tau( xtest, p-0, delta=0.8, hparm=2 )



args(ksm_tau)






# Checking code equivalence





















# START: hstar() --------------------------------
pn     <- hstar(abdord, wtord, const, n, tn)
tauscr <- ((asc[n] - asc[1]) * pn)/(2 * tn)
tauscr <- sqrt(n/(n - p)) * (1/tauscr)
w      <- rep(0, n)
stan   <- (ehat - median(ehat))/mad(ehat)
w[abs(stan) < hparm] <- 1
hubcor <- sum(w)/n
if (hubcor < 1e-06) {
  hubcor <- 1e-06
}
fincor <- 1 + (((p)/n) * ((1 - hubcor)/hubcor))
tauscr <- fincor * tauscr
tauscr



## START: ???????????????() -------------------




# Naive Attempt 1 --------------------------------

## Part 1: phi(1) - phi(0)
ordered_scores <- score_2s( 1:n, phi_2s_wilcox )
ordered_scores_deriv <- score_2s( 1:n, phi_2s_wilcox, dphi_2s_wilcox )
const <- ordered_scores[n] - ordered_scores[1]

eord <- sort(ehat)

## Part 2: c_n

osd <- ordered_scores_deriv
cn0 <- sum(osd)


tauss <- 0
for( ii in 1:n ){
  for( jj in 1:n ){
    it <- ordered_scores_deriv[jj] * ifelse( abs(eord[ii]-eord[jj]) <= 0.8, 1, 0 )
    tauss < -tauss + it
  }
}



outer( eord, eord, FUN="-")




# Check this constant thing ----------------------------------




ordered_scores <- score_2s( 1:n, phi_2s_wilcox )
ordered_scores_deriv <- score_2s( 1:n, phi_2s_wilcox, dphi_2s_wilcox )
const0 <- ordered_scores[n] - ordered_scores[1]

vphi0 <- ordered_scores_deriv
cn0 <- sum(ordered_scores_deriv)/const0


sum(ordered_scores_deriv)





pair_idx <- combn( n,2 )
vad0 <- cbind( ehat[pair_idx[1,]], ehat[pair_idx[2,]])
vwt0 <- cbind( vphi0[pair_idx[1,]], vphi0[pair_idx[2,]])/const0
vad02 <- abs(vad0[,1] - vad0[,2])
vwt02 <- vwt0[,1] + vwt0[,2]
iord <- order(vad02)
abs_diff_ord <- vad02[iord]
wt_ord <- vwt02[iord]



cumsum(wt_ord) / ((n - 1) * cn0)


## 

sum( 1:n )
sum(ordered_scores_deriv)

sum(wt_ord)



2*dphi/const0

















