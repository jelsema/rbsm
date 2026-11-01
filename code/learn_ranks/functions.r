
#'
#' Attempting to build my own R-estimators
#' 


# Replicate scores() function call -----------------

#' Will want some sort of selector function to set things like
#'   phi, dphi
#' and then return a function that produces both the scores and the 
#' derivative of the scores. Or maybe produces a function for each?



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
  theta <- optim( 
    par=c("theta"=median(x) ), fn=signedrank_dispersion, 
    method="Brent", lower=min(x), upper=max(x),
    sfun=sfun, pfun=pfun, x=x# ,
    # control=list( maxit = 200, reltol = .Machine$double.eps^(3/4) )
  )$par
  
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


# Regression -----------------

#' Page 168-169

jaeckel_dispersion <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  sum( dd*ss )
}

jaeckel_dispersion_deriv <- function( beta, x, y, sfun=score_2s, pfun=phi_2s_wilcox ){
  dd <- y - as.matrix(x) %*% beta
  ss <- sfun( rank(dd), pfun )
  ### sum( x * -1*ss )
  ### drop( t(x) %*% (-1*ss) )
  drop( crossprod( x, (-1*ss) ) )
}


phi_2s_wilcox <- function(u){ sqrt(12)*(u-0.5) }
dphi_2s_wilcox <- function(u){ rep(sqrt(12),length(u)) }
phi_2s_sign <- function(u){ sign(2*u -1) }
score_2s <- function( r, pfun, dpfun=NULL ){ 
  a0 <- pfun( r / (length(r)+1) )
  a1 <- a0 - mean(a0)
  a1sd <- sum(a1^2)/( length(a1)+1)
  scrs <- a1 / sqrt(a1sd)
  if( is.null(dpfun) ){
    return(scrs)
  } else{
    b0 <- dpfun( r / (length(r)+1) )
    return( b0/sqrt(a1sd) )
  }
}


# jdf_optim <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox ){
#   theta <- optim( par=c("theta"=lm(y~x)$coef[2] ), fn=jaeckel_dispersion, 
#                   method="Brent", lower=min(x), upper=max(x),
#                   sfun=sfun, pfun=pfun, y=y, x=x )$par
#   return( theta )
# }

jdf_optim <- function(y, x, sfun=score_2s, pfun=phi_2s_wilcox, method="Brent" ){
  
  # theta <- optim( par=c("theta"=lm(y~x)$coef[-1] ), fn=jaeckel_dispersion, 
  #                 method=method, lower=min(x), upper=max(x),
  #                 sfun=sfun, pfun=pfun, y=y, x=x )$par
  theta <- optim( par=c("theta"=lm(y~x)$coef[-1] ), fn=jaeckel_dispersion, 
                  method=method,
                  sfun=sfun, pfun=pfun, y=y, x=x )$par
  
  return( theta )
}


jaeckle_coef <- function( x, y, bhat0=NULL, yhat0=NULL, control=NULL, ... ){
  # Should add ability to pass score function
  # Calls to optim also need to have that added
  # Add symmetric or not as argument and to srd_optim
  
  X0 <- model.matrix( y ~ x )
  X1 <- X0[,-1,drop=FALSE]
  
  qrx <- qr(X0)
  Q0 <- qr.Q(qrx)
  XQ <- Q0[ , -1 , drop=FALSE ]
  
  # Initialize fit
  if( is.null( bhat0 ) ){
    if( is.null(yhat0) ){
      yhat0 <- y
    }
    bhat0 <- drop(crossprod( XQ , yhat0))
  }
  
  
  disp_null <- jaeckel_dispersion( 0*bhat0, XQ, y )
  disp_init <- jaeckel_dispersion( bhat0, XQ, y )
  if( disp_null < disp_init ){ bhat0 <- 0*bhat0 }
  
  # Optimize the jaeckle dispersion in the no-intercept space
  bhat_noint <- jaeckle_coef_noint( XQ, y, bhat0 )$par
  
  
  yhat <- XQ %*% bhat_noint
  ehat <- y - yhat
  
  alphahat <- srd_optim( ehat )
  ehat <- ehat - alphahat
  yhat <- yhat + alphahat
  # Should be able to replace lsfit() with a matrix / QR-approach
  # if something like lsfit doesn't exist in other languages
  bhat <- lsfit(X0, yhat, intercept = FALSE)$coef 
  # bhat
  names(bhat) <- NULL
  
  alphahat0 <- srd_optim(y)
  bhat_null <- c(alphahat0, rep(0, length(bhat) - 1))
  
  ehat1 <- y - X0 %*% bhat
  # tau_est <- gettauF0( ehat1, p=ncol(X1), scores = wscores )
  tau_est <- ksm_tau( ehat1, p=ncol(X1), ... )
  
  return(
    list(
      coefficients = bhat,
      residuals = ehat, fitted.values = yhat, 
      tau_est = tau_est
    )
  )
  
}

# cbind(
#   Rfit::rfit( Y ~ X, symmetric = TRUE )$coef,
#   jaeckle_coef( X, Y )$coef
# )



jaeckle_coef_noint <- function( x, y, bhat, control=NULL, ... ){
  # Follows the same logic/approach as Rfit
  # Ensure control has tolerance and iterations
  if( is.null(control) ){
    control <- list( maxit = 200, reltol = .Machine$double.eps^(3/4) )
  } else {
    if (is.null(control$reltol)) {
      control$reltol <- .Machine$double.eps^(3/4)
    }
    if (is.null(control$maxit)) {
      control$maxit <- 200
    }
  }
  
  ysd   <- sd(y)
  ys <- y/ysd
  
  fit0a <- optim(bhat/ysd, jaeckel_dispersion, method = "BFGS", 
                 x = x, y = ys,  gr = jaeckel_dispersion_deriv,
                 control = control
  )
  
  fit1a <- optim(fit0a$par * ysd, jaeckel_dispersion, method = "BFGS", 
                 x = x, y = y,  gr = jaeckel_dispersion_deriv,
                 control = control
  )
  
  return( fit1a )
}




# Tau -----------------------


ksm_tau <- function( ehat, p=0, delta=0.8, hparm=2, pfun=phi_2s_wilcox, dpfun=dphi_2s_wilcox ){
  
  n <- length(ehat)
  ordered_scores <- score_2s( 1:n, pfun )
  ordered_scores_deriv <- score_2s( 1:n, pfun, dpfun )
  
  # evaluate n vs p
  # if n/p >= 5, delta  = 0.80 works well
  # if n/p <  5, delta ~= 0.90 works well
  # Give warning for small n/p
  
  # START: hstarready / hstar_prepare
  # const0 == denominator in HM 3.7.1, converting phi to phi-star
  
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
  # END: hstarready / hstar_prepare
  
  # return( cbind(abs_diff_ord,wt_ord))  
  
  wt_ord_cum_nrmlz <- cumsum(wt_ord) / ((n - 1) * cn0)
  tn_out <- looptau0( abs_diff_ord, wt_ord_cum_nrmlz, delta )
  tn <- tn_out$tndelta / sqrt(n)
  pn <- hstar0(abs_diff_ord, wt_ord_cum_nrmlz, tn)
  # return( c(tn, pn) )
  
  gamma_est <- (const0 * pn)/(2 * tn)
  tau_est <- sqrt(n/(n-p)) * (1/gamma_est) # why not n-p-1 ?
  w <- rep(0, n)
  ehat_stnd <- (ehat - median(ehat))/mad(ehat)
  
  w[abs(ehat_stnd) < hparm] <- 1
  hubcor <- max( c(sum(w)/n, 1e-06) )
  # hubcor <- sum(w)/n
  # if (hubcor < 1e-06) {
  #   hubcor = 1e-06
  # }
  fincor = 1 + (((p)/n) * ((1 - hubcor)/hubcor))
  tau_est = fincor * tau_est
  return( tau_est )
}



# Replaces hstar()
hstar0 <- function( abs_diff_ord, wt_ord_cum_nrmlz, val ){
  # hstar_all <- cumsum(wt_ord) / ((n - 1) * cn)
  # Could perhaps use optim() here, but that'd likely slow the evaluation
  if( all( abs_diff_ord > val ) ){
    hstar_val <- 0
  } else{
    hstar_val <- max(wt_ord_cum_nrmlz[ abs_diff_ord <= val ])
  }
  return(hstar_val)
}


# Replaces looptau()
looptau0 <- function( abs_diff_ord, wt_ord_cum_nrmlz, delta ){
  
  m <- length(abs_diff_ord)
  idx_delta <- round(delta*m) # replace: inddelta
  tn_guess <- abs_diff_ord[idx_delta]  # replace: guess
  hval <- hstar0( abs_diff_ord, wt_ord_cum_nrmlz, tn_guess ) # replace: y
  
  max_iter <- m         # replace: itmax
  num_iter <- 0         # replace: ?icmax
  idx_keep <- idx_delta # replace: ikeep
  finished <- FALSE     # replace: gohome
  iError   <- 0
  
  if( hval == delta ){
    tn0 <- tn_guess
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
        tn0      <- tn_guess
        finished <- TRUE
        iError   <- 1
      } else{
        tn_guess = abs_diff_ord[idx_keep]
        hval = hstar0(abs_diff_ord, wt_ord_cum_nrmlz, tn_guess)
      }
    } else{
      # End loop if:
      # - moving left but hval < delta
      # - moving right but hval > delta
      tn0      <- tn_guess
      finished <- TRUE
      iError   <- 1
    }
    
    if( num_iter > max_iter ){
      tn0      <- tn_guess
      finished <- TRUE
      iError   <- 2
    }
  }
  
  return(
    list(
      tndelta = tn0,
      hval = hval, 
      iError = iError, 
      num_iter = num_iter, 
      finished = finished
    )
  )
  
}




