#' @title No-arbitrage panels
#' @description Calculate minimal distance no-arbitrage prices
#' @param option.panel A data-frame with columns relMid = normalized mid-prices (current stock-price = 1), r = annualized interest, logF = log forward price, dT = option maturity
#' @export
#' @return A no-arbitrage enforced panel with the same structure as the original panel
noArbPrices <- function(option.panel) {
  # get strikes
  k <- exp(option.panel$k)

  # turn everything into put prices
  put.price <- option.panel$relMid
  
  F <- exp(option.panel$logF[1])
  r <- option.panel$r[1]
  t <- option.panel$dT[1]
  q <- (-log(F)/t+r)
  
  put.price[k>F] <- put.price[k>F] - exp(-q*t) + exp(-r*t)*k[k>F]
  
  # set up quadratic programming, using the fact that the put prices are the integral of the distribution function
  x <- matrix(0,nrow = length(put.price),ncol=length(put.price))
  
  # the constant part
  x[,1] <- 1
  
  # the integral part
  for (nn1 in 1:(length(put.price)-1)) {
    for (nn2 in 1:nn1) {
      x[1+nn1,1+nn2] <- k[nn1+1]-k[nn2]
    }
  }
  
  # enforce non-negativity constraints on the INCREMENTS of the average distribution function values
  R <- diag(length(put.price))
  
  rvec <- rep(0,length(put.price))
  
  # run the L1 fitting of the options, given the no-arbitrage constraints
  no.arb <- rq.fit.fnc(x,put.price,R=R,r=rvec)
  
  fitted.prices <- x %*% no.arb$coefficients
  
  # calculate otm no-arbitrage prices
  otm.noarb <- fitted.prices
  otm.noarb[k>F] <- otm.noarb[k>F] + exp(-q*t) - exp(-r*t)*k[k>F]
  
  # replace mid-prices with no-arbitrage prices
  option.panel$relMid <- otm.noarb
  
  return(option.panel)
}