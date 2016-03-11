#' Calculate interpolated constant maturity implied GFTs
#' 
#' This function takes a list of priceMats and a markets, calculates the impleid GFT-s and then interpolatese/extrapolates them to the required maturities.
#' @param option.panels A list of Kx2 strike-price matrices
#' @param mkt.frame A data-frame showing the available maturities (corresponding to the day where the panels are observed) and corresponding r and q
#' @param u.seed A seed which is used to scale the frequencies by \code{\link{uSquareRootMat}}
#' @param mat.vec A vector of maturities where the interpolated gft-s are required.
#' @param offset A scalar/vector that specifies if we want to deviate from the square-root frequency calculation rule
#' @param time.IV A logical. If TRUE, then implied volatilities are time extrapolated from the thin-plate-spline. Otherwise transforms are calculated for maturities available in the option panels and the transforms are directly extrapolated.
#' @param ... Additional arguments passed to \code{\link{impliedGFT}}
#' @export
#' @return Returns an UxT matrix containing the interpolated gft at the scaled frequencies. If in u.seed we request to calculate portfolios for a u such that the power does not exist for too many portfolios (or some other error occurs in the \code{\link{spline}()} call), the value \code{-998} will be returned.

iGFTinterpolate <- function(option.panels, mkt.frame, u.seed, mat.vec, sell.offset=0, time.IV = TRUE, k.shrink = 1,...) {
  stopifnot(length(option.panels) == nrow(mkt.frame))
  
  T <- length(option.panels)
  
  # initialize outputs
  int.mat <- matrix(NA,length(u.seed), length(mat.vec))
  dimnames(int.mat)[[1]] <- paste("u",seq(1:length(u.seed)),sep=".")
  dimnames(int.mat)[[2]] <- paste("t",seq(1:length(mat.vec)),sep=".")
  int.mat.sell <- int.mat
  
  # we can either directly time extrapolate/interpolate from the tps, or we can do it on the transforms
  if (time.IV) {
    joined.vec <- sort(unique(c(mat.vec,mat.vec - sell.offset)))
    
    u.t.mat <- expand.grid(u = u.seed, t = joined.vec)
    
    # now calculate transforms
    transf.mat <- tpsImpliedGFT(option.panels,mkt.frame,u.t.mat,...)
      
    # now we use the linear time moment approximation for accurate interpolation
    for (uu in 1:length(u.seed)) {
      for (mm in 1:length(mat.vec)) {
        int.mat[uu,mm] <- transf.mat[transf.mat[,1] == u.seed[uu] & transf.mat[,2] == mat.vec[mm],3]
        int.mat.sell[uu,mm] <- transf.mat[transf.mat[,1] == u.seed[uu] & transf.mat[,2] == (mat.vec[mm] - sell.offset),3]
      }
    }
  
  } else {
    u.t.mat <- expand.grid(u = u.seed, t = mkt.frame$t)
    
    # now calculate transforms
    transf.mat <- tpsImpliedGFT(option.panels,mkt.frame,u.t.mat,...)
    
    # now we use the linear time moment approximation for accurate interpolation
    for (uu in 1:length(u.seed)) {
        int.mat[uu,] <- 1 + spline(mkt.frame$t,y=(Re(transf.mat[transf.mat[,1]==u.seed[uu],3])-1)/mkt.frame$t,xout=mat.vec,method="natural")$y * mat.vec
        int.mat.sell[uu,] <- 1 + spline(mkt.frame$t,y=(Re(transf.mat[transf.mat[,1]==u.seed[uu],3])-1)/mkt.frame$t,xout=(mat.vec - sell.offset),method="natural")$y * (mat.vec - sell.offset)
    }
    
  }
  
  return(list(gft.buy = int.mat, gft.sell = int.mat.sell))
}