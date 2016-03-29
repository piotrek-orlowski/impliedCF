#' Calculate the implied divergences from option panels
#' @param option.panels A list, each element containing a dataframe with normalized mid OTM option prices and LOG strikes. Each panel should correspond to the SAME day.
#' @param mkt.list A data-frame, where each row corresponds to the market features associated with the n-th option panel (e.g. maturity, interest rate, etc...)
#' @param u.t.mat An Nx3 dataframe, consisting of the powers, maturities, and divergence types
# #' @param doPlot Whether fitted / true values should be plotted
## ' @param doFitPlot Whether backfitting iterations should provide plot, \code{doFitPlot} is an integer: if 0, not plots, otherwise plots every \code{doFitPlot} iterations
#' @param spline.rank This is the \code{k} argument to the function \code{\link{s}} from package \link{mgcv}, it determines the spline basis dimension.
#' @param convergence.scale Scaling of the convergence criterion for the backfitting algorithm, a submodel \code{j} is deemed as having converged if the sqrt-sample-size-scaled norm of the iteration change in the predicted values is smaller than \code{convergence.scale * sd(cp - sum(k != j)fk)}
#' @param verbose Do we want diagnostic information?
#' @param ... further arguments to be passed to \code{\link{panelSmoother}} and further to \code{\link{rgl::persp3d}}
#' @export
#' @return An Nx3 dataframe, where the last column is the calculated transform value

tpsImpliedDivergence <- function(option.panels,mkt.frame,u.t.mat, doPlot=FALSE, doFitPlot = 0, verbose=FALSE, L = NULL, U = NULL, bootstrap = F, Nrepl = 1e3, ...) {
  
  # Fit model (with k=1)
  smoothing.results <- panelSmoother_asymptotic(option.panels, mkt.frame, ...) 
  tps.fit <- smoothing.results$model
  regr.mat <- smoothing.results$regr.mat
  otmFun <- cmpfun(smoothing.results$otmFun)
  otmFunCov <- cmpfun(smoothing.results$otmFunCov)
  
  # create output matrix
  u.t.out <- cbind(u.t.mat,res=NA)
  
  # now create a data-frame for interpolating interest-rates and dividend yields
  int.mat <- matrix(NA,length(unique(regr.mat$dT)),3)
  
  for (nn in 1:nrow(int.mat)) {
    int.mat[nn,1] <- unique(regr.mat$dT)[nn]
    int.mat[nn,2] <- subset(regr.mat,dT==int.mat[nn,1])$r[1]
    int.mat[nn,3] <- subset(regr.mat,dT==int.mat[nn,1])$q[1]
  }
  
  # create a matrix for storing prices simulated from the posterior parameter distribution
  sim.mat <- matrix(0,nrow(u.t.mat),1+Nrepl)
  sim.mat.quad <- sim.mat
  sim.par <- mvtnorm::rmvnorm(n = Nrepl, mean = coef(tps.fit), sigma = vcov(tps.fit))
  sim.par <- cbind(coef(tps.fit), t(sim.par))

  # quadrature nodes
  quad.nodes <- gauss.quad(n = 2^9, kind = "legendre")
  
  L0 <- L
  U0 <- U
  
  # now loop through maturity, frequency combinations and calculate implied transform
  for (nn in 1:nrow(u.t.mat)) {
    # calculate logF, q, and r
    t <- u.t.mat[nn,2]
    u <- u.t.mat[nn,1]
    
    r <- if(nrow(int.mat)>1){approx(int.mat[,1],int.mat[,2],xout=t,rule=2)$y} else {int.mat[,2]}
    q <- if(nrow(int.mat)>1){approx(int.mat[,1],int.mat[,3],xout=t,rule=2)$y} else {int.mat[,3]}
    
    if(inherits(smoothing.results$model,"gam")){
      s <- predict(smoothing.results$model, newdata = data.frame(k = (r-q)*t, dT = t)) 
    } else {
      s <- smoothing.results$model((r-q)*t)
    }
    
    if(is.null(L0)){
      L <- -10 * sqrt(t) 
    } else {
      L <- L0 * sqrt(t) 
    }
    if(is.null(U0)){
      U <- 10 * sqrt(t) 
    } else {
      U <- U0 * sqrt(t)
    }
    
    quad.nodes.loc <- quad.nodes
    quad.nodes.loc$nodes <- 0.5*(U-L)*quad.nodes.loc$nodes + 0.5*(U+L)
    
    logF <- (r-q)*t
    
    if(u.t.mat[nn,3] == 'div'){
      try({
        if(bootstrap){
          # sim.mat[nn,] <- adaptIntegrate(function(k) otmFunCov(tps.fit,k,r,q,t,sim.par) * matrix(exp(k*(u-1)) * exp(r*t), nrow=length(k), ncol = ncol(sim.par)), lowerLimit = L, upperLimit = U, fDim = ncol(sim.par), maxEval = 1e3)$integral
          otm.quad <- 0.5*(U-L)*otmFunCov(tps.fit,quad.nodes.loc$nodes,r,q,t,sim.par) * matrix(exp(quad.nodes.loc$nodes*(u-1)) * exp(r*t), nrow=length(quad.nodes.loc$nodes), ncol = ncol(sim.par))
          otm.quad <- otm.quad * matrix(quad.nodes.loc$weights, nrow = length(quad.nodes.loc$weights), ncol = ncol(sim.par))
          sim.mat[nn,] <- apply(otm.quad,2,sum)
        } else {
          u.t.out[nn,4] <- integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * exp(r*t), lower = logF, upper = U, subdivisions = 1000L, rel.tol = 1e-6)$value
          u.t.out[nn,4] <- u.t.out[nn,4] + integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * exp(r*t), lower = L, upper = logF, subdivisions = 1000L, rel.tol = 1e-6)$value
          # u.t.out[nn,4] <- u.t.out[nn,4] * exp(- u * (r-q)*t) 
        }
      })
    } else if(u.t.mat[nn,3] == 'skew'){
      try({
        if(bootstrap){
          # sim.mat[nn,] <- adaptIntegrate(function(k) otmFunCov(tps.fit,k,r,q,t,sim.par) * matrix(exp(k*(u-1)) * k * exp(r*t), nrow=length(k), ncol = ncol(sim.par)), lowerLimit = L, upperLimit = U, fDim = ncol(sim.par), maxEval = 1e3)$integral
          otm.quad <- 0.5*(U-L)*otmFunCov(tps.fit,quad.nodes.loc$nodes,r,q,t,sim.par) * matrix(exp(quad.nodes.loc$nodes*(u-1)) * quad.nodes.loc$nodes * exp(r*t), nrow=length(quad.nodes.loc$nodes), ncol = ncol(sim.par))
          otm.quad <- otm.quad * matrix(quad.nodes.loc$weights, nrow = length(quad.nodes.loc$weights), ncol = ncol(sim.par))
          sim.mat[nn,] <- apply(otm.quad,2,sum)
        } else {
          u.t.out[nn,4] <- integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * k * exp(r*t), lower = logF, upper = U, subdivisions = 1000L, rel.tol = 1e-6)$value
          u.t.out[nn,4] <- u.t.out[nn,4] + integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * k * exp(r*t), lower = L, upper = logF, subdivisions = 1000L, rel.tol = 1e-6)$value
          # u.t.out[nn,4] <- u.t.out[nn,4] * exp(- u * (r-q)*t) 
        }
      })
    } else if(u.t.mat[nn,3] == 'quart'){
      try({
        if(bootstrap){
          # sim.mat[nn,] <- adaptIntegrate(function(k) otmFunCov(tps.fit,k,r,q,t,sim.par) * matrix(exp(k*(u-1)) * k^2 * exp(r*t), nrow=length(k), ncol = ncol(sim.par)), lowerLimit = L, upperLimit = U, fDim = ncol(sim.par), maxEval = 1e3)$integral
          otm.quad <- 0.5*(U-L)*otmFunCov(tps.fit,quad.nodes.loc$nodes,r,q,t,sim.par) * matrix(exp(quad.nodes.loc$nodes*(u-1)) * quad.nodes.loc$nodes^2 * exp(r*t), nrow=length(quad.nodes.loc$nodes), ncol = ncol(sim.par))
      otm.quad <- otm.quad * matrix(quad.nodes.loc$weights, nrow = length(quad.nodes.loc$weights), ncol = ncol(sim.par))
      sim.mat[nn,] <- apply(otm.quad,2,sum)
        } else {
          u.t.out[nn,4] <- integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * k^2 * exp(r*t), lower = logF, upper = U, subdivisions = 1000L, rel.tol = 1e-6)$value
          u.t.out[nn,4] <- u.t.out[nn,4] + integrate(function(k) otmFun(tps.fit,k=k,r=r,q=q,dT=t) * exp(k*(u-1)) * k^2 * exp(r*t), lower = L, upper = logF, subdivisions = 1000L, rel.tol = 1e-6)$value
          # u.t.out[nn,4] <- u.t.out[nn,4] * exp(- u * (r-q)*t) 
        }
      })
    }
  }
  
  if(bootstrap){
    u.t.out[,4] <- apply(sim.mat,1,mean)
    return(list(iDiv = u.t.out, bSample = t(sim.mat)))
  } else {
    return(u.t.out)
  }
}