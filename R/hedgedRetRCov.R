#' Realized variance - covariance matrix of hedged option returns.
#' 
#' @description The function calculates the approximate matrix of realized variance-covariance between hedged GLT portfolios. The hedging frequency is determined by the frequency of stock futures data. If there is missing data in the time series of GLT values, these days are ommitted from hedging. The approximate nature of the result comes from the fact that we use the stochastic integral representation of the hedged GLT portfolio return and use Ito's isometry to build the approximation.
#' 
#' @note TODO: handle missing data (NA observations in the igft); add paper citation in @param psd.ensure
#' @param igft.list A list of buy and sell values of the igft (the difference comes from the maturity decreasing over the trading interval + vol. changes) as output by \code{\link{constMatiGFT}}
#' @param igft.days A vector containing the days corresponding to igft.list
#' @param hf.dframe A data frame containing high frequency futures prices. It has 2+#TTM columns, \code{day}, \code{time} and \code{F_ttm(1)...F_ttm(T)}, \code{T = length(fixed.maturities)}.
#' @param weekly.grid The approximately weekly grid of days where we empirically observe option prices (or their transformed values)
#' @param trade.end.grid (optional) A grid of days on which trades starting at weekly.grid are supposed to end. Use if weekly.grid has some 'gaps'.
#' @param u.seed Vector of base frequencies at which the portfolio values have been calculated. The vector's length is \code{N}
#' @param fixed.maturities The maturities of the option portfolios we are mimicking
#' @param psd.ensure If set to TRUE (default), the estimated variance covariance matrix is projected onto the space of positive semi-definite matrices as in [paper].
#' @param trade.int The trading interval (annualized in business year) over which we hold the option position.
#' @export
#' @return Array of size \code{(NxT)x(NxT)xS}, where \code{N} is the number of base frequencies, \code{T} is the number of maturities, \code{S} is the length of \code{weekly.grid}.
#' 

hedgedRetRCov <- function(igft.list, igft.days, hf.dframe, weekly.grid, trade.end.grid = NULL, u.seed, fixed.maturities = c(1/12,6/12), psd.ensure = TRUE, trade.int = 5/252){
  
  # check if the weekly grid is aligned with the option observations
  stopifnot(all(union(igft.days,weekly.grid) == igft.days))
  
  # we have one less return than observations
  weekly.ret.rCov <- array(NA,dim = c(length(u.seed)*length(fixed.maturities),length(u.seed)*length(fixed.maturities),length(weekly.grid)-1))
  
  # create effective u matrix from seed
  u.mat <- matrix(NA,length(u.seed),length(fixed.maturities))
  for (t in fixed.maturities) {u.mat[,which(t==fixed.maturities)] <- uSquareRootMat(u.seed,t)}
  
  print(u.mat)
  
  # now we create a function that calculates the appropriate realized variances and covariances
  weekly.rCov <- function(day.start,day.end,hf.stock.data) {
    # initialize array containing the hedged returns' vars and covars
    hedged.rCov <- array(NA,dim = c(length(u.seed)*length(fixed.maturities),length(u.seed)*length(fixed.maturities)))
    
    # get P0, the initial reference return
    P0 <- igft.list[[which(day.start==igft.days)]]$gft.buy
    P0[which(P0==-998)] <- NA
    
    # It might be easier to work if we don't vectorise in ttms, but let's take a look at that.
    for(ttm.1 in fixed.maturities){
      ttm.1.pos <- which(ttm.1 == fixed.maturities)
      for(ttm.2 in fixed.maturities){
        ttm.2.pos <- which(ttm.2 == fixed.maturities)
        for(uu.1 in u.seed){
          uu.1.pos <- which(uu.1 == u.seed)
          u.seed.less <- u.seed[which(uu.1 == u.seed):length(u.seed)]
          for(uu.2 in u.seed){
            uu.2.pos <- which(uu.2 == u.seed)
            # Prepare indexing
            row.index <- ttm.1.pos + uu.1.pos - 1
            col.index <- ttm.2.pos + uu.2.pos - 1
            
            # Prepare column names for selecting futures hedging data.
#             hedge.column.1 <- paste0("F_",sprintf("%1.4f",ttm.1))
#             hedge.column.2 <- paste0("F_",sprintf("%1.4f",ttm.2))
            # Temporary:
            hedge.column.1 <- hedge.column.2 <- 3
            
            # Check indexing -- you only want to do the lower triangular part of the variance - covariance matrix.
            if(col.index <= row.index){
              
              # For every pair (uu.1,uu.2) we want to calculate products of certain quantities. We first want to go through all days available for hedging and interpolate P_t,s intray-daily. Then put everything into vectors.
              
              options.hedgedays <- igft.days[igft.days >= day.start & igft.days <= day.end]
              # Initialise variables for holding anything you need. These will be vectors as long as your return hedge series
              hedge.vec.len <- nrow(subset(hf.stock.data,day %in% options.hedgedays[-1])) + 1
              P.tau.uu.1 <- P.tau.uu.2 <- rep(NA,hedge.vec.len - 1)
              
              for(dd in options.hedgedays[-1]){
                
#                 igft.buy <- igft.list[[which(dd == igft.days)]]$gft.buy
#                 igft.sell <- igft.list[[which(dd == igft.days)]]$gft.sell
                # Here we determine the time axis for interpolating the GLT values over the course of the day. The values in the first two columns of the hf.stock.data data.frame determine your approach to overnight changes and business time flow.
                
                # This determines the futures prices data frame for hedging. We assume that futures data are available for each day.
                hf.daily <- subset(hf.stock.data,day == dd)$F
                hf.daily.df <- subset(hf.stock.data,day == dd)
                
                # we are missing the overnight stock movement, so let's add that too
                hf.daily <- c(tail(subset(hf.stock.data,day<dd),1)$F,hf.daily)
                hf.daily.df <- rbind(tail(subset(hf.stock.data,day<dd),1),hf.daily.df)
                loc.time <- rowSums(hf.daily.df[,c(1,2),with = FALSE])
                
                dd.prev <- tail(options.hedgedays[which(options.hedgedays < dd)],1)
                if(dd.prev == day.start){
                  igft.prev <- igft.list[[which(dd.prev == igft.days)]]$gft.buy
                }
                else{
                  igft.prev <- igft.list[[which(dd.prev == igft.days)]]$gft.sell 
                }
                
                igft.now <- igft.list[[which(dd == igft.days)]]$gft.sell
                
                # This function works for igft.buy and igft.sell to be values at the start and end of the day. It assumes we have observations for days in igft.days, which come more often than observing returns via weekly.grid. 
                
                dframe.index <- which(hf.stock.data$day == dd)
                dframe.index <- c(head(dframe.index,1)-1 , dframe.index[-length(dframe.index)])
                dframe.index <- dframe.index - (nrow(subset(hf.stock.data,day == options.hedgedays[which(options.hedgedays == dd)-1])) - 1)
                
                P.tau.uu.1[dframe.index] <- approx(x = c(loc.time[1],loc.time[length(loc.time)]), y = Re(c(igft.prev[uu.1.pos,ttm.1.pos],igft.now[uu.1.pos,ttm.1.pos])), xout = loc.time[-length(loc.time)])$y + 1i * approx(x = c(loc.time[1],loc.time[length(loc.time)]), y = Im(c(igft.prev[uu.1.pos,ttm.1.pos],igft.now[uu.1.pos,ttm.1.pos])), xout = loc.time[-length(loc.time)])$y
                
                if((uu.2 == uu.1) & (ttm.1 == ttm.2)){
                  P.tau.uu.2[dframe.index] <- P.tau.uu.1[dframe.index]
                }
                else{
                  P.tau.uu.2[dframe.index] <- approx(x = c(loc.time[1],loc.time[length(loc.time)]), y = Re(c(igft.prev[uu.2.pos,ttm.2.pos],igft.now[uu.2.pos,ttm.2.pos])), xout = loc.time[-length(loc.time)])$y + 1i * approx(x = c(loc.time[1],loc.time[length(loc.time)]), y = Im(c(igft.prev[uu.2.pos,ttm.2.pos],igft.now[uu.2.pos,ttm.2.pos])), xout = loc.time[-length(loc.time)])$y  
                }
                
                # you have to multiply each P.tau by the appropriate power of the stock price.
                
                P.tau.uu.1[dframe.index] <- P.tau.uu.1[dframe.index] * hf.stock.data[dframe.index,hedge.column.1]^uu.1
                if(uu.2 == uu.1){
                  P.tau.uu.2[dframe.index] <- P.tau.uu.1[dframe.index]
                }
                else{
                  P.tau.uu.2[dframe.index] <- P.tau.uu.2[dframe.index] * hf.stock.data[dframe.index,hedge.column.2]^uu.1
                }  
              }
            }
            # Now you have long time series of P and hedging terms. You can directly compute the cross products.
            temp.cov <- P.tau.uu.1[-length(P.tau.uu.1)]*P.tau.uu.2[-length(P.tau.uu.2)]/(P.tau.uu.1[1]*P.tau.uu.2[1])
            temp.cov <- temp.cov * (diff(P.tau.uu.1)/P.tau.uu.1[-length(P.tau.uu.1)] - uu.1 * diff(hf.stock.data[,hedge.column.1])/hf.stock.data[-nrow(hf.stock.data),hedge.column.1])
            temp.cov <- temp.cov * (diff(P.tau.uu.2)/P.tau.uu.2[-length(P.tau.uu.2)] - uu.2 * diff(hf.stock.data[,hedge.column.2])/hf.stock.data[-nrow(hf.stock.data),hedge.column.2])
            hedged.rCov[row.index, col.index] <- sum(temp.cov)
          }
        }  
      }  
    }
    #print(paste("Hedged ret:",hedged.ret))
    
    # Produce a Hermitian matrix
    hedged.rCov <- hedged.rCov + t(Conj(hedged.rCov)) - diag(diag(Re(hedged.rCov)))
    
    # Some diagnostics of the matrix.
    # First, hedged.rCov might contain complex values off the diagonal.
    # Calculate the eigenvalues to check if it is positive definite:
    rCov.eig <- eigen(hedged.rCov, symmetric = TRUE)
    if(min(rCov.eig$values) <0){
      warning(paste0("Non-positive definite variance-covariance matrix for returns starting at ",day.start))
    }
    if((min(rCov.eig$values) < 0 ) & psd.ensure){
      vec <- t(rCov.eig$vectors)
      val <- rCov.eig$values
      val[which(val < 0)] <- 0
      hedged.rCov <- t(vec) %*% diag(val) %*% vec
    }
    
    return(hedged.rCov)
  }
  
  # convert hf data-frame to data table for faster subsetting
  hf.tab <- data.table(hf.dframe, key="day")
  
  # if the weekly.grid is indeed close to weekly, do nothing. If it has gaps, use trade.end.grid to determine day.end.
  if(is.null(trade.end.grid)){
    for (ww in 1:(length(weekly.grid)-1)) {
      if (ww %% 300 == 0) {
        print(ww)
      }
      # load appropriate hf futures returns
      hf.stock.data <- hf.tab[day >= weekly.grid[ww] & day <= weekly.grid[ww+1],]
      weekly.ret.cov[,,ww] <- weekly.rCov(weekly.grid[ww],weekly.grid[ww+1], hf.stock.data)
      if (ww %% 300 == 0) {
#         print(apply(weekly.ret[,,1:ww,drop=FALSE],c(1,2),mean,na.rm=TRUE))
#         print(apply(weekly.ret[,,1:ww,drop=FALSE],c(1,2),var,na.rm=TRUE))
      }
    }  
  }
  else{
    for (ww in 1:(length(weekly.grid)-1)) {
      if (ww %% 300 == 0) {
        print(ww)
      }
      # load appropriate hf futures returns
      hf.stock.data <- hf.tab[day >= weekly.grid[ww] & day <= trade.end.grid[ww],]
      weekly.ret.rCov[,,ww] <- weekly.rCov(weekly.grid[ww],trade.end.grid[ww], hf.stock.data)
      if (ww %% 300 == 0) {
#         print(apply(weekly.ret.rCov[,,1:ww,drop=FALSE],c(1,2),mean,na.rm=TRUE))
#         print(apply(weekly.ret.rCov[,,1:ww,drop=FALSE],c(1,2),var,na.rm=TRUE))
      }
    }
  }
  
  # check if you have to return complex or real
  if(max(abs(Im(weekly.ret.rCov)),na.rm=TRUE)<=1e-5){
    weekly.ret.rCov <- Re(weekly.ret.rCov)
  }
  
  return(list(weekstart = weekly.grid[-length(weekly.grid)],hedged.cov = weekly.ret.cov))
}