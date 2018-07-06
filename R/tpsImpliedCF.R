#' Calculate the implied characteristic function using a joint fitting of a smooth surface to all options.
#' @param option.panels A list, each element containing a dataframe with normalized mid OTM option prices and LOG strikes. Each panel should correspond to the SAME day.
#' @param mkt.list A data-frame, where each row corresponds to the market features associated with the n-th option panel (e.g. maturity, interest rate, etc...)
#' @param u.t.mat An Nx2 dataframe, consisting of the imaginary part of the argument and maturities for which the CFs should be calculated
#' @param discounted logical, if \code{TRUE}, the discounted CF (price of CF) is returned, otherwise the RN cf.
#' @param doPlot Whether fitted / true values should be plotted
#' @param doFitPlot Whether backfitting iterations should provide plot, \code{doFitPlot} is an integer: if 0, not plots, otherwise plots every \code{doFitPlot} iterations
#' @param spline.rank This is the \code{k} argument to the function \code{\link{s}} from package \link{mgcv}, it determines the spline basis dimension.
#' @param convergence.scale Scaling of the convergence criterion for the backfitting algorithm, a submodel \code{j} is deemed as having converged if the sqrt-sample-size-scaled norm of the iteration change in the predicted values is smaller than \code{convergence.scale * sd(cp - sum(k != j)fk)}
#' @param weights string which defines the weighting scheme to be used in fitting. \code{opt-sqrt} weighs by the square root of the out of the money option price.
#' @param verbose Do we want diagnostic information?
#' @param ... further arguments to be passed to \code{\link{panelSmoother}} and further to \code{\link{rgl::persp3d}}
#' @export
#' @return An Nx3 dataframe, where the last column is the calculated transform value
tpsImpliedCF <- function(option.panels, mkt.frame, u.t.mat, discounted = TRUE, smoothing.results = NULL, doPlot=FALSE, doFitPlot = 0, verbose=FALSE, ...) {
  
  # Fit model (with k=1)
  if(is.null(smoothing.results)){
    smoothing.results <- panelSmoother(option.panels, mkt.frame, ...) 
  }
  tps.fit <- smoothing.results$model
  regr.mat <- smoothing.results$regr.mat
  otmFun <- cmpfun(smoothing.results$otmFun)
  
  if (doPlot>0) {
    t <- unique(regr.mat$dT)[doPlot]
    pan.sub <- subset(regr.mat,dT==t)
    k.pred <- seq(12*quantile(pan.sub$k,1e-4),12*quantile(pan.sub$k,1-1e-4),length.out=1000)
    otm.pred <- otmFun(tps.fit,k=k.pred,r=pan.sub$r[1],q=pan.sub$q[1],dT=t)
    plot(k.pred,predict(tps.fit,data.frame(k=k.pred,dT=t)),type="l")
    points(pan.sub$k,pan.sub$IV,col="red")
  }
  
  # create output matrix
  u.t.out <- cbind(u.t.mat,NA)
  
  # now create a data-frame for interpolating interest-rates and dividend yields
  int.mat <- matrix(NA,length(unique(regr.mat$dT)),3)
  
  for (nn in 1:nrow(int.mat)) {
    int.mat[nn,1] <- unique(regr.mat$dT)[nn]
    int.mat[nn,2] <- subset(regr.mat,dT==int.mat[nn,1])$r[1]
    int.mat[nn,3] <- subset(regr.mat,dT==int.mat[nn,1])$q[1]
  }
  
  # now loop through maturity, frequency combinations and calculate implied transform
  for (nn in 1:nrow(u.t.mat)) {
    # calculate logF, q, and r
    t <- u.t.mat[nn,2]
    u <- u.t.mat[nn,1]
    
    r <- approx(int.mat[,1],int.mat[,2],xout=t,rule=2)$y
    q <- approx(int.mat[,1],int.mat[,3],xout=t,rule=2)$y
    
    s <- predict(smoothing.results$model, newdata = data.frame(k = 0, dT = t))
    
    L <- -10 * sqrt(t) * s
    U <- 10 * sqrt(t) * s
    
    logF <- (r-q)*t
    
  #   try({u.t.out[nn,3] <- cos(u*(r-q)*t) + 1i*sin(u*(r-q)*t);
  #   u.t.out[nn,3] <- u.t.out[nn,3] + 
  #     integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*cos(u*k)+u*sin(u*k))/exp(k) * exp(r*t),lower=logF,upper=U,subdivisions = 20000L, rel.tol = 1e-8)$value + 
  #     1i*integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*sin(u*k)-u*cos(u*k))/exp(k) * exp(r*t),lower=logF,upper=U,subdivisions = 20000L, rel.tol = 1e-8)$value
  #   u.t.out[nn,3] <- u.t.out[nn,3] + 
  #     integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*cos(u*k)+u*sin(u*k))/exp(k) * exp(r*t),lower=L,upper=logF,subdivisions = 20000L, rel.tol = 1e-8)$value + 
  #     1i*integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*sin(u*k)-u*cos(u*k))/exp(k) * exp(r*t),lower=L,upper=logF,subdivisions = 20000L, rel.tol = 1e-8)$value})
  #   u.t.out[nn,3] <- u.t.out[nn,3] * exp(- u*(r-q)*t*1i)
  # }
    try({u.t.out[nn,3] <- cos(u*(r-q)*t) + 1i*sin(u*(r-q)*t);
      u.t.out[nn,3] <- u.t.out[nn,3] + 
        pcubature(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*cos(u*k)+u*sin(u*k))/exp(k) * exp(r*t),lowerLimit = logF,upperLimit = U, tol = 1e-6)$integral + 
        1i*pcubature(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*sin(u*k)-u*cos(u*k))/exp(k) * exp(r*t),lowerLimit = logF,upperLimit = U,tol = 1e-6)$integral
      u.t.out[nn,3] <- u.t.out[nn,3] + 
        pcubature(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*cos(u*k)+u*sin(u*k))/exp(k) * exp(r*t),lowerLimit = L,upperLimit = logF, tol = 1e-6)$integral + 
        1i*pcubature(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * (-u^2*sin(u*k)-u*cos(u*k))/exp(k) * exp(r*t),lowerLimit = L,upperLimit = logF, tol = 1e-6)$integral})
      u.t.out[nn,3] <- u.t.out[nn,3] * ifelse(discounted,exp(- u*(r-q)*t*1i),1)
  }
  return(u.t.out)
}