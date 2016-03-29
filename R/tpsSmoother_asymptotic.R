#' @export
panelSmoother_asymptotic <- function(option.panels, mkt.frame, gam.bs = 'ds', gam.m = c(1,0.5), ...){
  
  # function for rescaling IV for asymtptotic behaviour
  # asFoo <- function(x,t,tr) {x <- x/sqrt(t);(abs(x) <= tr) *(  log(1+tr)  ) +(abs(x) > tr) * log(1+abs(x))}
  # asFoo <- function(x,t) {x <- x/sqrt(t);(abs(x) <= 0.3) *(  log(1+0.3)  ) +(abs(x) > 0.3) * log(1+abs(x))}
  asFoo <- function(x,t) log(2+abs(x))
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
    
    put <- call <- o.pan.noarb$relMid
    
    put[o.pan[,"k"]>logF] <- -exp(-q*dT) + exp(o.pan[,"k"][o.pan[,"k"]>logF] - r*dT) + put[o.pan[,"k"]>logF]
    call[o.pan[,"k"]<=logF] <- exp(-q*dT) - exp(o.pan[,"k"][o.pan[,"k"]<=logF] - r*dT) + call[o.pan[,"k"]<=logF]
    
    IV <- sapply(1:length(call),function(x) tryCatch(vanillaOptionImpliedVol(exercise="european",price=call[x],S=1,X=exp(o.pan[,"k"][x]),tau=dT,r=r,q=q), error=function(e) NA_real_))
    # w <- vanillaOptionEuropean(S=1,X=exp(o.pan[,"k"]),tau=dT,r=r,q=q,v=IV)$vega
    w <- sapply(1:length(call), function(x) tryCatch(vanillaOptionEuropean(S=1,X=exp(o.pan[x,"k"]),tau=dT,r=r,q=q,v=IV[x])$vega, error=function(e) NA_real_))
    
    w <- w^2/sum(w^2)
    
    # save time-to-maturity, moneyness and log put+call
    kvec <- o.pan[,"k"]-logF
    rownames(regr.mat) <- NULL
    regr.mat <- rbind(regr.mat,cbind(dT,kvec,IV,w,r,q))
    regr.mat <- regr.mat[which(!is.na(regr.mat[,"IV"])),]
  }
  rownames(regr.mat) <- NULL
  regr.mat <- as.data.frame(regr.mat)
  colnames(regr.mat) <- c("dT","k","IV","w","r","q")
  
  tr <- max(abs(regr.mat$k/sqrt(regr.mat$dT))) + 0.05
  
  # regr.mat <- cbind(regr.mat, scIV = regr.mat$IV/asFoo(regr.mat$k,regr.mat$dT,tr))
  regr.mat <- cbind(regr.mat, scIV = regr.mat$IV/asFoo(regr.mat$k,regr.mat$dT))
  
  if(length(unique(regr.mat$dT)) == 1){
    tps.fit <- tryCatch(
      gam(IV ~ s(k,m=1,bs=gam.bs),data=regr.mat)
      , error = function(e){
        res <- tryCatch(gam(IV ~ s(k,m=1,bs="tp"),data=regr.mat),error = function(e){approxfun(x = regr.mat$k, y = regr.mat$IV, rule = 2)})
        return(res)
      }
    )
  } else {
    # tps.fit <- gam(IV ~ s(k,dT,m=c(1,0.5),bs="ds"),data=regr.mat)
    tps.fit <- gam(scIV ~ s(k,dT,m=gam.m,bs=gam.bs),data=regr.mat)  
  }
  
  
  otmFun <- function(models,k,r,q,dT){
    logF <- (r-q)*dT
    kvec <- (k-logF)
    
    # Create prediction data matrix to be used with the models -- remember that they were fit to de-meaned data.
    pred.mat <- data.frame(k=kvec,dT=dT)
    if(inherits(models,"gam")){
      IV.pred <- pmax(1e-4,predict(models, newdata = pred.mat))
      # IV.pred <- IV.pred * asFoo(x = kvec, t = dT, tr = max(abs(regr.mat$k/sqrt(regr.mat$dT)))+0.05)
      IV.pred <- IV.pred * asFoo(x = kvec, t = dT)
      IV.pred <- pmax(IV.pred, 1e-2)
      IV.pred <- IV.pred^2
    } else {
      IV.pred <- pmax(1e-4, models(pred.mat[,"k"]))^2
    }
    
    put <- vanillaOptionEuropean(S=1,X=exp(k[kvec<=0]),tau=dT,r=r,q=q,v=IV.pred[kvec<=0],greeks=FALSE,type="put")
    call <- vanillaOptionEuropean(S=1,X=exp(k[kvec>0]),tau=dT,r=r,q=q,v=IV.pred[kvec>0],greeks=FALSE,type="call")
    return(c(put,call))
  }
  
  otmFunCov <- function(models,k,r,q,dT,coeffVec){
    logF <- (r-q)*dT
    kvec <- (k-logF)
    
    # Create prediction data matrix to be used with the models -- remember that they were fit to de-meaned data.
    pred.mat <- data.frame(k=kvec,dT=dT)
    if(inherits(models,"gam")){
      IV.pred <- predict(object = models, newdata = pred.mat, type = "lpmatrix")
      IV.pred <- IV.pred %*% coeffVec
      # IV.pred <- IV.pred * matrix(asFoo(x = kvec, t = dT, tr = max(abs(regr.mat$k/sqrt(regr.mat$dT)))+0.05), nrow = length(kvec), ncol = ncol(IV.pred), byrow = F)
      IV.pred <- IV.pred * matrix(asFoo(x = kvec, t = dT), nrow = length(kvec), ncol = ncol(IV.pred), byrow = F)
      IV.pred <- pmax(IV.pred, 1e-2)
      IV.pred <- IV.pred^2
    } else {
      IV.pred <- pmax(1e-4, models(pred.mat[,"k"]))^2
    }
    
    put <- apply(IV.pred,2,function(v.vec) vanillaOptionEuropean(S=1,X=exp(k[kvec<=0]),tau=dT,r=r,q=q,v=v.vec[kvec<=0],greeks=FALSE,type="put"))
    call <- apply(IV.pred,2,function(v.vec) vanillaOptionEuropean(S=1,X=exp(k[kvec>0]),tau=dT,r=r,q=q,v=v.vec[kvec>0],greeks=FALSE,type="call"))
    
    return(rbind(put,call))
  }
  
  return(list(model = tps.fit, regr.mat = regr.mat, otmFun = otmFun, otmFunCov = otmFunCov))
  
}