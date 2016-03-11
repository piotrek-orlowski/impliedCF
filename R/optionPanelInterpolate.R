#' Fit a smooth IV surface to all available options and return option prices for given maturities and strike ranges that depend on maturity and IV.
#' @param option.panels A list, each element containing a dataframe with normalized mid OTM option prices and LOG strikes. Each panel should correspond to the SAME day.
#' @param mkt.list A data-frame, where each row corresponds to the market features associated with the n-th option panel (e.g. maturity, interest rate, etc...)
#' @param out.mat Sx2 matrix, where first column specifies target maturities and the second specifies the number of options available at each maturity (for now they have to be equal)
#' @param doPlot Whether fitted / true values should be plotted
#' @param verbose Do we want diagnostic information?
#' @export
#' @return out.panels list with fields \code{p.list} -- a list of option panels, \code{m.list} -- a list of corresponding market structures.
optionPanelInterpolate <- function(option.panels,mkt.frame,out.mat, doPlot=FALSE, doFitPlot = 0, verbose=FALSE) {
  
  smoothing.results <- panelSmoother(option.panels, mkt.frame)
  tps.fit <- smoothing.results$model
  regr.mat <- smoothing.results$regr.mat
  otmFun <- smoothing.results$otmFun
  
  # create output objects
  p.list <- list()
  
    
  # Prepare mkt structure
  dT <- unique(regr.mat$dT)
  q <- vapply(dT,FUN.VALUE=numeric(1),FUN=function(x) return(head(subset(regr.mat,dT==x)$q,1)))
  r <- vapply(dT,FUN.VALUE=numeric(1),FUN=function(x) return(head(subset(regr.mat,dT==x)$r,1)))
  mkt <- data.frame(p = 1, q = approx(dT,q,xout=out.mat[,2],rule= 2)$y, r = approx(dT,r,xout=out.mat[,2],rule= 2)$y, t = out.mat[,2])
  
  if(any(is.na(mkt$r))){
    browser()
  }
  
  # Prepare strike grids
  
  strike.bounds.data <- vapply(option.panels,FUN.VALUE=numeric(2), FUN = function(x) return(range(x[,"k"])))
  
  # At some maturities, it seems there only is the call or put side, for some reason. Work around that by removing these data points from strike.bounds.data
  narrow.str <- which(strike.bounds.data[1,] >=0 & strike.bounds.data[2,] <= 0)
  if(length(narrow.str) > 0){
    strike.bounds.data <- strike.bounds.data[,-narrow.str]
  }
  
  # Interpolate strike grids to required maturities
  
  strike.bounds.projected <- -exp(spline(mkt.frame$t,log(-strike.bounds.data[1,]),xout = mkt$t)$y)
  strike.bounds.projected <- rbind(strike.bounds.projected,exp(spline(mkt.frame$t,log(strike.bounds.data[2,]),xout = mkt$t)$y))
  
  # In case the strike grids are all identical...
  
  strike.bounds.projected <- strike.bounds.projected + matrix(runif(prod(dim(strike.bounds.projected)), -1e-6,1e-6),nrow(strike.bounds.projected),ncol(strike.bounds.projected))
  
  # Create strike grids for the interpolated prices
  
  strike.grids <- apply(strike.bounds.projected, 2, function(x){
    out <- seq(x[1],x[2],length.out = out.mat[which(x[1] == strike.bounds.projected[1,]),1])
    return(out) 
  })
  
  # Calculate option prices
  option.prices <- apply(strike.grids, 2, function(x){
    return(otmFun(models = tps.fit, k = x, r = mkt$r[which(x[1] == strike.bounds.projected[1,])], q = mkt$q[which(x[1] == strike.bounds.projected[1,])], dT = mkt$t[which(x[1] == strike.bounds.projected[1,])]))
  })
  
  for(kk in 1:ncol(option.prices)){
    p.list[[kk]] <- matrix(c(strike.grids[,kk],option.prices[,kk]),ncol=2)
    colnames(p.list[[kk]]) <- c("k","relMid")
  }
  
  return(list(p.list = p.list, mkt = mkt))
}