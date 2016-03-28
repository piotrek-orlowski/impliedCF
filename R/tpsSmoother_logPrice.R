#' @export
panelSmoother_logPrice <- function(option.panels, mkt.frame, convex = F, gam.bs = 'ds', gam.m = c(1,0.5), ...){
  
  # create regression frame from list of option panels
  regr.mat <- NULL
  opt.mat <- NULL
  for (nn in 1:length(option.panels)) {
    o.pan <- option.panels[[nn]]
    
    opt.mat <- rbind(opt.mat, matrix(c(o.pan[,"k",drop=F],o.pan[,"relMid",drop=F]),ncol=2))
    
    if (median(o.pan[,"k",drop=F]) > .5) {warning("Seems as if surfaceImpliedGFT was supplied with effective strikes instead of LOG strikes...!")}
    
    # retrieve market features
    mkt <- mkt.frame[nn,]
    q <- mkt$q
    r <- mkt$r
    dT <- mkt$t
    logF <- (r-q) * dT
    
    # apply no-arbitrage filter
    if(nrow(o.pan) > 3){
      o.pan.noarb <- noArbPrices(option.panel=data.frame(k=o.pan[,"k"], relMid=o.pan[,"relMid"],r=r,dT=dT,logF=logF))
    } else {
      o.pan.noarb <- data.frame(k = o.pan[,"k"], relMid = o.pan[,"relMid"], r = r, dT = dT, logF = logF)
    }
    
    put <- call <- o.pan.noarb[,"relMid"]
    
    put[o.pan[,"k"]>logF] <- -exp(-q*dT) + exp(o.pan[,"k"][o.pan[,"k"]>logF] - r*dT) + put[o.pan[,"k"]>logF]
    call[o.pan[,"k"]<=logF] <- exp(-q*dT) - exp(o.pan[,"k"][o.pan[,"k"]<=logF] - r*dT) + call[o.pan[,"k"]<=logF]
    
    IV <- sapply(1:length(call),function(x) tryCatch(vanillaOptionImpliedVol(exercise="european",price=call[x],S=1,X=exp(o.pan[,"k"][x]),tau=dT,r=r,q=q), error=function(e) NA_real_))
    # w <- vanillaOptionEuropean(S=1,X=exp(o.pan[,"k"]),tau=dT,r=r,q=q,v=IV)$vega
    w <- sapply(1:length(call), function(x) tryCatch(vanillaOptionEuropean(S=1,X=exp(o.pan[x,"k"]),tau=dT,r=r,q=q,v=IV[x])$vega, error=function(e) NA_real_))
    
    w <- w^2/sum(w^2)
    
    o.pan.noarb[,"relMid"] <- log(o.pan.noarb[,"relMid"])
    # save time-to-maturity, moneyness and log put+call
    # kvec <- ((o.pan[,"k"]-logF)/sqrt(dT)^1)
    kvec <- o.pan[,"k"]-logF
    rownames(regr.mat) <- NULL
    regr.mat <- rbind(regr.mat,cbind(dT,kvec,IV,w,r,q,relMid = o.pan.noarb[,"relMid"], logPut = log(put), logCall = log(call)))
    regr.mat <- regr.mat[which(!is.na(regr.mat[,"IV"])),]
  }
  
  rownames(regr.mat) <- NULL
  regr.mat <- as.data.frame(regr.mat)
  colnames(regr.mat) <- c("dT","k","IV","w","r","q","relMid","logPut","logCall")
  
  if(length(unique(regr.mat$dT)) == 1){
    tps.fit <- tryCatch(
      gam(IV ~ s(k,m=1,bs=gam.bs),data=regr.mat)
      , error = function(e){
        res <- tryCatch(gam(IV ~ s(k,m=1,bs="tp"),data=regr.mat),error = function(e){approxfun(x = regr.mat$k, y = regr.mat$IV, rule = 2)})
        return(res)
      }
    )
  } else {
    if(convex){
      mod.scam.put <- scam(formula = logPut ~ s(k,bs="micv",m=2) + s(dT,bs="micv",m=2), data = regr.mat)
      mod.scam.call <- scam(formula = logCall ~ s(k,bs="mdcv",m=2) + s(dT,bs="micv",m=2), data = regr.mat)
    } else {
      mod.scam.put <- scam(formula = logPut ~ s(k,bs="mpi",m=2) + s(dT,bs="mpi",m=2), data = regr.mat)
      mod.scam.call <- scam(formula = logCall ~ s(k,bs="mpd",m=2) + s(dT,bs="mpi",m=2), data = regr.mat) 
    }
  }
  
  # otmFun for model that predicts log prices of calls and puts
  otmFun <- function(models,k,r,q,dT){
    newdata <- data.frame(k=k,r=r,q=q,dT=dT)
    pred.put <- predict(models$modput, newdata)
    pred.call <- predict(models$modcall, newdata)
    pred.out <- pred.call
    pred.out[which(newdata$k <= (newdata$r-newdata$q)*newdata$dT)] <- pred.put[which(newdata$k <= (newdata$r-newdata$q)*newdata$dT)]
    pred.out <- exp(pred.out)
    return(pred.out)
  }
  
  # otmFun for model that predicts log prices of calls and puts
  otmFunCov <- function(models,k,r,q,dT,coeffVec){
    newdata = data.frame(k=k,r=r,q=q,dT=dT)
    pred.call <- predict(models$modput, newdata, type = 'lpmatrix')
    pred.put <- predict(models$modcall, newdata, type = 'lpmatrix')
    pred.call <- pred.call %*% coeffVec
    pred.put <- pred.put %*% coeffVec
    pred.out <- pred.call
    pred.out[which(newdata$k <= (newdata$r-newdata$q)*newdata$dT)] <- pred.put[which(newdata$k <= (newdata$r-newdata$q)*newdata$dT)]
    pred.out <- exp(pred.out)
    return(pred.out)
  }

  return(list(model = list(modput = mod.scam.put, modcall = mod.scam.call), regr.mat = regr.mat, otmFun = otmFun, otmFunCov = otmFunCov))
  
}