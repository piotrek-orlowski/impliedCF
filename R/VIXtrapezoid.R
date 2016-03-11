#' VIX calculation using a trapezoid quadrature
#' 
#' This function calculates the value of the VIX using a simple trapezoid integral (and a Black-Scholes correction)
#' @param priceMat The Nx2 strike + otm price matrix
#' @param mkt The usual market list
#' @param cutoff.strike length-2 vector, determines the options to be discarded with strikes too low (\code{<= cutoff.strike[1]}) or too high (\code{>= cutoff.strike[2]})
#' @param cutoff.price numeric, options with price below this value will be discarded
#' @export
#' @return Returns the estimated value of the VIX for the given maturity
VIXtrapezoid <- function(priceMat,mkt,cutoff.strike = c(min(priceMat[,1]),max(priceMat[,1])), cutoff.price = 1e-7) {
  require('NMOF')
  
  stopifnot(cutoff.strike[2] > cutoff.strike[1])
  stopifnot(cutoff.strike[2] > 1)
  stopifnot(cutoff.strike[1] < 1)
  
  # remove prices with strikes too far out and those that are too low to bother
  if(is.null(colnames(priceMat))){
    colnames(priceMat) <- c("k","p")
  }
  priceMat <- subset(priceMat, priceMat[,"k"] >= cutoff.strike[1] & priceMat[,"k"] <= cutoff.strike[2])
  priceMat <- subset(priceMat, priceMat[,"p"] >= cutoff.price)
  
  # check if you're left with a sufficient number of options
  stopifnot(nrow(priceMat) >= 10)
  if(nrow(priceMat) <= 15){
    warning("You only have 15 options in your panel, possibly due to using cutoff.strike and cutoff.price, beware of results")
  }
  
  kvec <- priceMat[,1]
  
  # make sure strikes are ordered
  stopifnot(all(diff(kvec)>0))
  price.corrected <- price.otm <- priceMat[,2]
  
  # calculate futures value
  F <- exp(mkt$t*(mkt$r-mkt$q))
  
  # calculate ATM vol
  last.put.ind <- sum(kvec <= F)
  first.call.ind <- last.put.ind + 1
  
  IVbelow <- vanillaOptionImpliedVol(exercise="european",price=price.otm[last.put.ind],S=1,X=kvec[last.put.ind],tau=mkt$t,r=0,q=0,type="put")
  IVabove <- vanillaOptionImpliedVol(exercise="european",price=price.otm[first.call.ind],S=1,X=kvec[first.call.ind],tau=mkt$t,r=0,q=0,type="call")
  
  # we will use this volatility to correct the prices (the best would be the ATM vol, but we aim for simplicity here)
  vol.corr <- approx(kvec[c(last.put.ind,first.call.ind)],c(IVbelow,IVabove),F)$y * sqrt(mkt$t)
  
  # substract Black-Scholes prices to avoid the "ATM peak" problem
  price.corrected[kvec<=F] <- price.otm[kvec<=F] - vanillaOptionEuropean(S=F,X=kvec[kvec<=F],tau=mkt$t,r=0,q=0,v=vol.corr^2/mkt$t,greeks=FALSE,type="put")
  price.corrected[kvec>F] <- price.otm[kvec>F] - vanillaOptionEuropean(S=F,X=kvec[kvec>F],tau=mkt$t,r=0,q=0,v=vol.corr^2/mkt$t,greeks=FALSE,type="call")
  
  # apply trapezoid approximation to the corrected prices
  vix <- 2*sum((price.corrected[-1]/kvec[-1]^2 + price.corrected[-length(price.corrected)]/kvec[-length(kvec)]^2)/2 * diff(kvec)) * exp(mkt$r*mkt$t)
  
  vix <- vix + vol.corr^2
  return(vix/mkt$t)
}