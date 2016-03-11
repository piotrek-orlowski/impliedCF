#' Calculate the implied GFT for a series of option panels, using a square-root frequency scaling and for fixed maturities.
#' 
#' This function takes as input a vector of strings corresponding to option panels and then calculates the implied GFT for a given frequency rule. Then the implied GFT is interpolated for fixed maturities.
#' @param panels.list A list version of the option panels, as outputted by \code{\link{panelsToList}}
#' @param u.seed A vector of frequencies that will be scaled according to a square-root time-to-maturity rule
#' @param fixed.maturities A vector of maturities which represent the option trading strategies
#' @param trade.int A scalar representing the assumed trading time length. This together with the fixed.maturities defines the interpolating maturities.
#' @param output.name The file name where the interpolated igft-s will be saved
#' @param doParallel If integer, specifies the number of cores to be used by starting a cluster with the \code{parallel} package. If \code{foreach}, an extant parallel backend will be used with the \code{foreach} package.
#' @param ... further arguments to be passed to \code{\link{iGFTinterpolate}} and further on to \code{\link{tpsImpliedGFT}}, \code{\link{panelSmoother}} and further to \code{\link{rgl::persp3d}}
#' @export
constMatiGFT <- function(panels.list, u.seed, fixed.maturities = c(1/12,6/12), trade.int, output.name, doParallel = 2,...) {
  n.days <- length(panels.list$days.all)

  igft.list <- list()
  
  require('parallel')
  if(as.character(doParallel) == "foreach"){
    # Check for backend, stop if none
    backend.registered <- getDoParRegistered()
    stopifnot(backend.registered)
  }
  else{
    cl <- makeCluster(doParallel)
    
    clusterExport(cl,c("fixed.maturities","trade.int","u.seed"),envir=environment())  
    clusterEvalQ(cl,{library("impliedCF")}) 
  }
  
  # create joined list of panels and markets
  joined.list <- vector("list",length=(length(panels.list$p.list)))
  for (nn in 1:length(panels.list$p.list)) {joined.list[[nn]] <- list(pn = panels.list$p.list[[nn]], m =panels.list$m.list[[nn]])}
  
  calc.igft <- function(j.list) {
    i.pol <- iGFTinterpolate(j.list$pn, mkt.frame=j.list$m, u.seed, fixed.maturities, sell.offset = trade.int,...)
    print(i.pol)
    return(i.pol)
  }
  
  if(as.character(doParallel) == "foreach"){
    igft.list <- foreach(list = joined.list, .packages = 'impliedCF') %dopar% calc.igft(list)
  } else if (doParallel > 1) {
    igft.list <- parLapply(cl,joined.list,calc.igft)
  } else {
    igft.list <- lapply(joined.list,calc.igft)
  }
  if(as.character(doParallel) == "foreach"){
    # You don't care  
  }else{
    stopCluster(cl)
  }
  
  
  days.all <- panels.list$days.all
  save(igft.list, days.all, u.seed, fixed.maturities, trade.int, file=output.name)
  print(igft.list[[1]])
  return(0)
}