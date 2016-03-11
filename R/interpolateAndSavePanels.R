#' Calculate the implied GFT for a series of option panels, using a square-root frequency scaling and for fixed maturities.
#' 
#' This function takes as input a vector of strings corresponding to option panels and then calculates the implied GFT for a given frequency rule. Then the implied GFT is interpolated for fixed maturities.
#' @param panels.list A list version of the option panels, as outputted by \code{\link{panelsToList}}
#' @param out.mat Sx2 matrix, where first column specifies target maturities and the second specifies the number of options available at each maturity (for now they have to be equal).
#' @param output.name The file name where the interpolated igft-s will be saved
#' @param doParallel If integer, specifies the number of cores to be used by starting a cluster with the \code{parallel} package. If \code{foreach}, an extant parallel backend will be used with the \code{foreach} package.
#' @param k.shrink what?!
#' @export
saveInterpolatedPanels <- function(panels.list, out.mat, output.name, doParallel = 2, lib.dir = "test.lib.dir",...) {
  
  n.days <- length(panels.list$days.all)
  
  require('parallel')
  if(as.character(doParallel) == "foreach"){
    # Check for backend, stop if none
    backend.registered <- getDoParRegistered()
    stopifnot(backend.registered)
  }
  else{
    cl <- makeCluster(doParallel)
    
    clusterExport(cl,c("lib.dir","out.mat"),envir=environment())  
    clusterEvalQ(cl,{library("impliedCF",lib.loc=lib.dir)}) 
  }
  
  # create joined list of panels and markets
  joined.list <- vector("list",length=(length(panels.list$p.list)))
  for (nn in 1:length(panels.list$p.list)) {joined.list[[nn]] <- list(pn = panels.list$p.list[[nn]], m =panels.list$m.list[[nn]])}
  
  doInterpolate <- function(j.list) {
    for(kk in 1:length(j.list$pn)){
      j.list$pn[[kk]] <- data.frame(k = j.list$pn[[kk]][,1],relMid = j.list$pn[[kk]][,2])
#       colnames(j.list$pn[[kk]]) <- c("k","relMid")
    }
    i.pol <- optionPanelInterpolate(j.list$pn, j.list$m, out.mat)
    return(i.pol)
  }
  if(as.character(doParallel) == "foreach"){
    fitted.opts <- foreach(list = joined.list, .packages = 'impliedCF') %dopar% doInterpolate(list)
  } else if (doParallel > 1) {
    fitted.opts <- parLapply(cl,joined.list,doInterpolate)
  } else {
    fitted.opts <- lapply(joined.list,doInterpolate)
  }
  if(as.character(doParallel) == "foreach"){
    # You don't care  
  }else{
    stopCluster(cl)
  }
  
  # Convert back to lists that we like -- p.list, m.list, days.all
  
  fitted.opts.ret <- list(p.list = list(), m.list = list(), days.all = panels.list$days.all)
  fitted.opts.ret$p.list <- lapply(fitted.opts, function(ll) {
    return(ll$p.list)
    })
  fitted.opts.ret$m.list <- lapply(fitted.opts, function(ll) {
    return(ll$mkt)
  })
  
  days.all <- panels.list$days.all
  save(fitted.opts.ret, days.all, out.mat, file=output.name)
  print(fitted.opts.ret[[1]][[1]])
  return(0)
}