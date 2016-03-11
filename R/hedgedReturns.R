#' Hedged option returns
#' 
#' This function calculates the hedged option trading returns, using high frequency option data and approximately daily calculated igft.
#' @description The function calculates returns to trading hedged GLT positions, given a grid of GLT price data and a time series of high frequency stock futures data. The hedging frequency is determined by the frequency of stock futures data. If there is missing data in the time series of GLT values, these days are ommitted from hedging.
#' @note TODO: Adapt the function so that old calls work as intended, but new calls use daily data for calculating hedged returns; instead of dropping missing-data days from hedging, do wider interpolation; Ensure that with multiple maturities, multiple time series of futures prices are used.
#' @param igft.list A list of buy and sell values of the igft (the difference comes from the maturity decreasing over the trading interval + vol. changes) as output by \code{\link{constMatiGFT}}
#' @param igft.days A vector containing the days corresponding to igft.list
#' @param hf.dframe A data frame containing high frequency futures returns. Has 3 columns, day, time and F (which is the futures price)
#' @param weekly.grid The approximately weekly grid of days where we empirically observe option prices (or their transformed values)
#' @param trade.end.grid (optional) A grid of days on which trades starting at weekly.grid are supposed to end. Use if weekly.grid has some 'gaps'.
#' @param fixed.maturities The maturities of the option portfolios we are mimicking
#' @param trade.int The trading interval (annualized in business year) over which we hold the option position.
#' @export
#' @return Returns a data frame containing the starting day of each option position (assumed to be bought at the end of the day) and the corresponding weekly hedged return. 
hedgedReturns <- function(igft.list, igft.days, hf.dframe, weekly.grid, trade.end.grid = NULL, u.seed, fixed.maturities = c(1/12,6/12), trade.int = 5/252) {
  # check if the weekly grid is aligned with the option observations
  stopifnot(all(union(igft.days,weekly.grid) == igft.days))
  
  # we have one less return than observations
  weekly.ret <- array(NA,c(length(weekly.grid)-1,length(u.seed),length(fixed.maturities)))
  
  # create effective u matrix from seed
  u.mat <- matrix(u.seed,length(u.seed),length(fixed.maturities))  
  
  print(u.mat)
  
  # now we create a function that calculates the appropriate weekly return
  weekly.hedged.ret <- function(day.start,day.end,hf.stock.data) {
    # initialize matrix containing the hedged returns
    hedged.ret <- matrix(NA,length(u.seed),length(fixed.maturities))
    
    # now add unhedged return. If the return is not available for a given freq and date (code -998), NAs will be propagated and an NA will be returned at the end.
    P0 <- igft.list[[which(day.start==igft.days)]]$gft.buy
    Pt <- igft.list[[which(day.end==igft.days)]]$gft.sell

    P0 <- P0 * head(hf.stock.data$F,1)^(u.mat)
    Pt <- Pt * tail(hf.stock.data$F,1)^(u.mat)
    hedged.ret <- Pt / P0  - 1

    # now let's look at the days where we are having our portfolio and for which option observations are available
    options.hedgedays <- igft.days[igft.days >= day.start & igft.days < day.end]
    
    # now iterate over the days and calculate the hedging portfolio values (need time interpolation). We assume option buying occurs at the END of the first day and selling occurs at the END of the last day.
    for (dd in options.hedgedays) {
      # let's get daily hf futures data. For this we need the current time interval (what days happened since last option observation)
      
      day.upper <- c(options.hedgedays[-1],day.end)[which(dd==options.hedgedays)]
      day.lower <- dd

      hf.daily <- subset(hf.stock.data,day>day.lower & day<=day.upper)$F
      
      # we are missing the overnight stock movement, so let's add that too
      hf.daily <- c(tail(subset(hf.stock.data,day<=day.lower),1)$F,hf.daily)

      igft.buy <- igft.list[[which(dd == igft.days)]]$gft.buy
      igft.sell <- igft.list[[which(dd == igft.days)]]$gft.sell

      # let's interpolate the igft to correspond to the actual option maturity left
      igft.ipolated <- matrix(NA,length(u.seed),length(fixed.maturities))
      for (uu in 1:length(u.seed)) {
        
        for (tt in 1:length(fixed.maturities)) {
          # note that everything should be in business time!
          igft.ipolated[uu,tt] <- approx(x=c(day.start,day.end),y=Re(c(igft.buy[uu,tt],igft.sell[uu,tt])),xout=dd)$y + 1i * approx(x=c(day.start,day.end),y=Im(c(igft.buy[uu,tt],igft.sell[uu,tt])),xout=dd)$y 
          
          P.ipolated <- igft.ipolated[uu,tt] * head(hf.daily,-1)^u.mat[uu,tt]
          delta.ipolated <- P.ipolated / head(hf.daily,-1)
          
          # now add return on each hedging trade
          hedged.ret[uu,tt] <- hedged.ret[uu,tt] - u.mat[uu,tt] * sum(delta.ipolated * diff(hf.daily)) / (P0[uu,tt])
        }
      }
    }
    return(hedged.ret)
  }
  
  # convert hf data-frame to data table for faster subsetting
  hf.tab <- data.table(hf.dframe, key="day")
  
  # if the weekly.grid is indeed close to weekly, do nothing. If it has gaps, use trade.end.grid to determine day.end.
  if(is.null(trade.end.grid)){
    for (ww in 1:(length(weekly.grid)-1)) {
      if (ww %% 1 == 0) {
        print(ww)
      }
      # load appropriate hf futures returns
      hf.stock.data <- hf.tab[hf.tab$day > weekly.grid[ww] & hf.tab$day <= weekly.grid[ww+1],]
      hf.stock.data <- rbind(tail(hf.tab[hf.tab$day <= weekly.grid[ww],],1),hf.stock.data)
      weekly.ret[ww,,] <- weekly.hedged.ret(weekly.grid[ww],weekly.grid[ww+1], hf.stock.data)
      if (ww %% 10 == 0) {
        print(apply(weekly.ret[1:ww,,,drop=FALSE],c(2,3),mean,na.rm=TRUE))
        print(apply(weekly.ret[1:ww,,,drop=FALSE],c(2,3),var,na.rm=TRUE))
      }
    }  
  }
  else{
    for (ww in 1:(length(weekly.grid)-1)) {
      if (ww %% 300 == 0) {
        print(ww)
      }
      # load appropriate hf futures returns
      hf.stock.data <- hf.tab[day > weekly.grid[ww] & day <= trade.end.grid[ww],]
      weekly.ret[ww,,] <- weekly.hedged.ret(weekly.grid[ww],trade.end.grid[ww], hf.stock.data)
      if (ww %% 300 == 0) {
        print(apply(weekly.ret[1:ww,,,drop=FALSE],c(2,3),mean,na.rm=TRUE))
        print(apply(weekly.ret[1:ww,,,drop=FALSE],c(2,3),var,na.rm=TRUE))
      }
    }
  }
  
  
  # check if you have to return complex or real
  if(max(abs(Im(weekly.ret)),na.rm=TRUE)<=1e-5){
    weekly.ret <- Re(weekly.ret)
  }
  
  return(list(weekstart = weekly.grid[-length(weekly.grid)],hedged.ret = weekly.ret))
}