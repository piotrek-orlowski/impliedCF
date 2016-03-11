#' Calculate the implied Laplace transform using a joint fitting of a smooth surface to all options.
#' @param option.panels A list, each element containing a dataframe with normalized mid OTM option prices and LOG strikes. Each panel should correspond to the SAME day.
#' @param mkt.list A data-frame, where each row corresponds to the market features associated with the n-th option panel (e.g. maturity, interest rate, etc...)
#' @param u.t.mat An Nx2 dataframe, consisting of the frequencies and maturities for which the laplace transforms should be calculated
#' @param doPlot Whether fitted / true values should be plotted
#' @param doFitPlot Whether backfitting iterations should provide plot, \code{doFitPlot} is an integer: if 0, not plots, otherwise plots every \code{doFitPlot} iterations
#' @param spline.rank This is the \code{k} argument to the function \code{\link{s}} from package \link{mgcv}, it determines the spline basis dimension.
#' @param convergence.scale Scaling of the convergence criterion for the backfitting algorithm, a submodel \code{j} is deemed as having converged if the sqrt-sample-size-scaled norm of the iteration change in the predicted values is smaller than \code{convergence.scale * sd(cp - sum(k != j)fk)}
#' @param weights string which defines the weighting scheme to be used in fitting. \code{opt-sqrt} weighs by the square root of the out of the money option price.
#' @param verbose Do we want diagnostic information?
#' @export
#' @return An Nx3 dataframe, where the last column is the calculated transform value
surfaceImpliedGFT <- function(option.panels,mkt.frame,u.t.mat, max.pow = 3, doPlot=FALSE, doFitPlot = 0, verbose=FALSE, spline.rank = 6, convergence.scale = 1e-2, weights = "opt-sqrt") {
  
  # create regression frame from list of option panels
  regr.mat <- NULL
  opt.mat <- NULL
  for (nn in 1:length(option.panels)) {
    o.pan <- option.panels[[nn]]
    
    opt.mat <- rbind(opt.mat, matrix(c(o.pan$k,o.pan$relMid),ncol=2))
    
    if (median(o.pan$k) > .5) {warning("Seems as if surfaceImpliedGFT was supplied with effective strikes instead of LOG strikes...!")}
      
    put <- call <- o.pan$relMid
    
    # retrieve market features
    mkt <- mkt.frame[nn,]
    q <- mkt$q
    r <- mkt$r
    dT <- mkt$t
    logF <- (r-q) * dT
    
    put[o.pan$k>logF] <- -exp(-q*dT) + exp(o.pan$k[o.pan$k>logF] - r*dT) + put[o.pan$k>logF]
    call[o.pan$k<=logF] <- exp(-q*dT) - exp(o.pan$k[o.pan$k<=logF] - r*dT) + call[o.pan$k<=logF]

    # save time-to-maturity, moneyness and log put+call
    regr.mat <- rbind(regr.mat,cbind(dT,(o.pan$k-logF)/sqrt(dT)^1,log(put)+log(call),r,q))    
  }
  regr.mat <- as.data.frame(regr.mat)
  
  colnames(regr.mat) <- c("dT","k","cp","r","q")
  colnames(opt.mat) <- c("k","otm")
  
  # weight each maturity "equally" // this is now overriden two lines below, constant weighting is used.
  weight.vec <- regr.mat$dT
  for (nn in 1:length(weight.vec)) {weight.vec[nn] <- 1/(sum(regr.mat$dT == weight.vec[nn]))}
  weight.vec <- rep(1,nrow(regr.mat))

  # You will be backfitting an additive model where the linear parametric part has linear inequality restrictions on the parameters.  
  
  # What variables will enter the prediction? Create formula string for the linear model.
  pred.vars <- c("abs(k)","k","1")
  pred.str <- NULL
  pred.indiv <- NULL
  pred.indiv.index <- 1
  t.powers <- c(0,1,2:max.pow)
  for (pp in pred.vars) {
    for (mm in t.powers) {
      pred.indiv[pred.indiv.index] <- paste0("I(",pp,"*dT^",mm,")")
      pred.indiv.index <- pred.indiv.index + 1
      pred.str <- paste(pred.str,paste0("I(",pp,"*dT^",mm,")"),sep="+")
    }
  }
  
  pred.indiv <- c("cp",pred.indiv)
  
  # Prepare formula for linear model fitting, removing the constant.
  formula.str <- paste0("cp ~ -1 ", pred.str)
  
  # De-mean data (cp only).
#   regr.mat.means <- apply(as.matrix(regr.mat), 2, mean)
  regr.mat.means <- rep(0,length(colnames(regr.mat)))
  regr.mat.means[which(colnames(regr.mat) == "cp")] <- mean(regr.mat[,"cp"])
  names(regr.mat.means) <- colnames(regr.mat)
  regr.mat.copy <- regr.mat
  regr.mat <- regr.mat - matrix(regr.mat.means, nrow(regr.mat), ncol(regr.mat), byrow = TRUE)
  
  # Enforce negative terms in k^2 and abs(k) for integration stability.
  t.vec <- unique(u.t.mat[,2])
  t.vec <- t.vec - regr.mat.means["dT"]
#   R <- matrix(0,length(t.vec)*2,length(pred.vars) * (max.pow + 1))
#   R <- matrix(0,length(t.vec)*3,length(pred.vars) * (max.pow + 1))
  R <- matrix(0,length(t.vec)*3,length(pred.vars) * (length(t.powers)))
#   R <- rbind(R,matrix(0,length(t.vec), length(pred.vars) * length(t.powers)))
#   R <- matrix(0,length(t.vec)*5,length(pred.vars) * (length(t.powers)))
  
  for (nn in 1:length(t.vec)) {
    for (mm in t.powers) {
      mm.index <- which(mm == t.powers)-1
      R[nn,1+mm.index] <- -t.vec[nn]^mm # absolute values
      R[nn+1*length(t.vec),c(1+mm.index, 1*length(t.powers)+1+mm.index)] <- rep(-t.vec[nn]^mm,2)# abs and k on the positive half-line
      R[nn+2*length(t.vec),c(1+mm.index, 1*length(t.powers)+1+mm.index)] <- c(-t.vec[nn]^mm, t.vec[nn]^mm) # abs and k on the negative half-line
    }
  }
#   for(nn in 1:length(t.vec)){
#     for(mm in which((t.powers-2)>=0)){
#       R[nn+4*length(t.vec), 3*length(t.powers) + mm] <- - t.powers[mm]*(t.powers[mm]-1)*t.vec[nn]^(max(0,t.powers[mm]-2))
#     }
#   }
  
  r <- 0 * R[,1] + 1e-6
#   r[(length(t.vec)+1):(2*length(t.vec))] <- -5
  r[(1*length(t.vec)+1):(2*length(t.vec))] <- abs(max(u.t.mat[,1]))
  r[(2*length(t.vec)+1):(3*length(t.vec))] <- abs(min(u.t.mat[,1]))

  # Fit initial linear model.
  meq <- 0
  option.lm.base <- constr.lm.estim(pred.indiv, regr.mat, R, r, weight.vec, meq = meq)
  
#   option.lm.base <- orlm(model = lm(as.formula(formula.str), data = regr.mat, weights = weight.vec), ui = R, ci = r, index = 1:ncol(R), orig.out = TRUE)
  
  # Introduce smoothing models weights, overridden, weights are now equail
#   gam.wt <- 1/pmax(abs(regr.mat$k + regr.mat.means["k"]),1e-2)^2*weight.vec^2
  gam.wt <- rep(1,nrow(regr.mat))
  
  # make a list of model formulae
  gam.form.list <- list()
  gam.form.list[[1]] <- NULL
  gam.form.list[[2]] <- "cp ~ s(dT,k = min(length(option.panels)-1))"
  gam.form.list[[3]] <- paste0("cp ~ s(k, m = 1, k = ", spline.rank,")")
  gam.form.list[[4]] <- paste0("cp ~ s(I(k*dT), m = 1, k = ", spline.rank,")")
  gam.form.list[[5]] <- paste0("cp ~ s(I(abs(k) * dT), m = 1, k = ", spline.rank,")")
  gam.form.list[[6]] <- paste0("cp ~ s(I(k^2*dT^(0.5)), m = 1, k = ", spline.rank,")")
  
  # Fit initial smoothing models, one by one, while removing the fitted values from already fitted models -- this is the 0th iteration of the backfitting algorithm.
  regr.mat.loc <- regr.mat
  regr.mat.loc$cp <- regr.mat.loc$cp - option.lm.base$fitted.values
  option.gam.1.base <- gam(cp ~ s(dT,k = min(length(option.panels)-1)), data = regr.mat.loc, weights = gam.wt)
  regr.mat.loc$cp <- regr.mat.loc$cp - option.gam.1.base$fitted.values
  option.gam.2.base <- gam(cp ~ s(k, m = 1, k = spline.rank), data = regr.mat.loc, weights = gam.wt)
  regr.mat.loc$cp <- regr.mat.loc$cp - option.gam.2.base$fitted.values
  option.gam.3.base <- gam(cp ~ s(I(k*dT), m = 1, k = spline.rank), data = regr.mat.loc, weights = gam.wt)
  regr.mat.loc$cp <- regr.mat.loc$cp - option.gam.3.base$fitted.values
  option.gam.4.base <- gam(cp ~ s(I(abs(k) * dT), m = 1, k = spline.rank), data = regr.mat.loc, weights = gam.wt)
  regr.mat.loc$cp <- regr.mat.loc$cp - option.gam.4.base$fitted.values
  option.gam.5.base <- gam(cp ~ s(I(k^2 * dT^0.5), m = 1, k = spline.rank), data = regr.mat.loc, weights = gam.wt)
  
  model.list.old <- list(option.lm.base, option.gam.1.base, option.gam.2.base, option.gam.3.base, option.gam.4.base, option.gam.5.base)
  
  # With the residuals from the initial fit, define weights
  if(weights == "cp-relative"){
    cp.fitted <- rowSums(matrix(vapply(model.list.old[1:5], FUN = function(x) return(x$fitted.values), FUN.VALUE = rep(0,nrow(regr.mat))),ncol=5))
    gam.wt <- weight.vec <- pmax(1/(pmax(abs(cp.fitted - regr.mat$cp)/abs(regr.mat$cp),0.15)), 1)^2
    weight.vec <- sqrt(weight.vec)
  }
  else if(weights == "cp-absolute"){
    cp.fitted <- rowSums(matrix(vapply(model.list.old[1:5], FUN = function(x) return(x$fitted.values), FUN.VALUE = rep(0,nrow(regr.mat))),ncol=5))
    gam.wt <- weight.vec <- pmax(1/pmax(abs(cp.fitted - regr.mat$cp),0.02),0.2)^2
    weight.vec <- sqrt(weight.vec)
  }
  else if(weights == "opt-sqrt"){
#     gam.wt <- weight.vec <- pmax(opt.mat[,"otm"]/sqrt(regr.mat[,"dT"]),1e-1)
#     weight.vec <- pmax(sqrt(opt.mat[,"otm"]/sqrt(regr.mat[,"dT"])), sqrt(1e-1))
    gam.wt <- weight.vec <- pmax(opt.mat[,"otm"],1e-3)
    weight.vec <- pmax(sqrt(opt.mat[,"otm"]), sqrt(1e-3))
#     
  }
  else{
    gam.wt <- weight.vec <- rep(1,nrow(regr.mat))
  }
  # Initialise variables necessary for backfitting
  convergence.fit <- FALSE
  .iter <- 0
  model.vec <- c(1,2,3,4,5,6)
#   model.vec <- c(1,2,3,4,5)
  unconverged.models <- model.vec
  model.list.new <- model.list.old

  # Initialise variables necessary for plotting in the fitting process.
  plotcnt <- 0
  colvec <- rainbow(9)

  # Start the fitting exercise
  while(!convergence.fit & .iter < 200){
    .iter <- .iter + 1
    for(kk in unconverged.models){
      # Remove all-except-one models' predictions from the data:
      regr.mat.loc <- regr.mat
      regr.mat.loc$cp <- regr.mat$cp - rowSums(vapply(model.list.new[setdiff(model.vec,kk)], FUN = function(x){return(x$fitted.values)}, FUN.VALUE = rep(0,nrow(regr.mat))))
      
      if(kk == 1){
#         tmp <- orlm(model = lm(as.formula(formula.str), data = regr.mat.loc, weights = weight.vec), ui = R, ci = r, index = 1:ncol(R), orig.out = TRUE)
        model.list.new[[kk]] <- constr.lm.estim(terms = pred.indiv, data = regr.mat.loc, R, r, weight.vec, meq=meq)
      }
      else{
        model.list.new[[kk]] <- gam(as.formula(gam.form.list[[kk]]), data = regr.mat.loc, weights = gam.wt)
        model.list.new[[kk]]$fitted.values <- model.list.new[[kk]]$fitted.values - mean(model.list.new[[kk]]$fitted.values)
      }
      
    }
    
    # This function checks for convergence and, in case of verbose=TRUE, reports on fitting progress
    convergence.fit <- apply(matrix(model.vec), 1, function(m){
      loc.norm <- base::norm(matrix(model.list.old[[m]]$fitted.values - model.list.new[[m]]$fitted.values), type="F") / sqrt(length(model.list.new[[m]]$fitted.values))
      if(verbose){
        print(paste("In iteration ",.iter,"in model number",m,"the discrepancy is",loc.norm)) 
      }
      return(loc.norm < max((convergence.scale * sd(regr.mat$cp - rowSums(vapply(model.list.new[setdiff(model.vec,m)], FUN = function(x){return(x$fitted.values)}, FUN.VALUE = rep(0,nrow(regr.mat)))))),.Machine$double.eps))
    })
    
    # Here you can generate data and fitted plots every doFitPlot iterations.
    if(doFitPlot > 0){
      if(.iter %% doFitPlot == 0){
        if(plotcnt == 0){
          plot(regr.mat$cp,col="black")
        }
        plotcnt <- plotcnt + 1
        points(rowSums(vapply(model.list.new[setdiff(model.vec,12)], FUN = function(x){return(x$fitted.values)}, FUN.VALUE = rep(0,nrow(regr.mat)))),col=colvec[plotcnt %% 9 + 1],pch=4)
      } 
    }
    
    # Check which models converged and stop iterating over them
#     unconverged.models <- model.vec[-which(convergence.fit)]
    unconverged.models <- model.vec

    # Check whether all models converged.
    convergence.fit <- all(convergence.fit)
    model.list.old <- model.list.new
  }
  
  if (verbose) {
    for(kk in 1:length(model.list.new)){
      print(summary(model.list.new[[kk]]))
    }
  }
  
  # create function for OTM option price calculation from fitted GAM
  otmFun <- function(models,k,r,q,dT){
    b <- -exp(k - r*dT) + exp(-q*dT)
    logF <- (r-q)*dT
    kvec <- (k-logF)/sqrt(dT)^1
#     
    # Create prediction data matrix to be used with the models -- remember that they were fit to de-meaned data.
    pred.mat <- data.frame(k=kvec - regr.mat.means["k"],dT=dT - regr.mat.means["dT"])
    
    cp.pred <- rowSums(matrix(vapply(models[2:length(model.vec)], FUN = predict, FUN.VALUE = rep(0,nrow(pred.mat)), newdata = pred.mat), ncol=(length(model.vec)-1)))
    cp.pred <- cp.pred + predict(models[[1]], newdata = pred.mat)
    cp.pred <- cp.pred + regr.mat.means["cp"]
    c <- -exp(cp.pred)
    
    # calculate the two roots
    r1 <- 1/2*(-b+sqrt(b^2-4*c))
    r2 <- 1/2*(-b-sqrt(b^2-4*c))
    put <- pmax(r1,r2)
    put[put < .Machine$double.eps] <- 0
    call <- put+b
    call[call < .Machine$double.eps] <- 0
    return(c(put[kvec<=0],call[kvec>0]))
  }
  # create output matrix
  u.t.out <- cbind(u.t.mat,NA)
  
  # now create a data-frame for interpolating interest-rates and dividend yields
  int.mat <- matrix(NA,length(unique(regr.mat$dT)),3)
  
  for (nn in 1:nrow(int.mat)) {
    int.mat[nn,1] <- unique(regr.mat.copy$dT)[nn]
    int.mat[nn,2] <- subset(regr.mat.copy,dT==int.mat[nn,1])$r[1]
    int.mat[nn,3] <- subset(regr.mat.copy,dT==int.mat[nn,1])$q[1]
  }

  # now loop through maturity, frequency combinations and calculate implied transform
  for (nn in 1:nrow(u.t.mat)) {
    # calculate logF, q, and r
    t <- u.t.mat[nn,2]
    u <- u.t.mat[nn,1]
    
    r <- approx(int.mat[,1],int.mat[,2],xout=t,rule=2)$y
    q <- approx(int.mat[,1],int.mat[,3],xout=t,rule=2)$y
    
    L <- -10 * sqrt(t)
    U <- 10 * sqrt(t)
    
    logF <- (r-q)*t
    try({u.t.out[nn,3] <- exp(u*(r-q)*t);
        u.t.out[nn,3] <- u.t.out[nn,3]+integrate(function(k) otmFun(model.list.new,k=k,r=r,q=q,dT=t) * u * (u-1) * exp((u-1)*k) * exp(r*t),lower=logF,upper=U,subdivisions = 100L)$value
        u.t.out[nn,3] <- u.t.out[nn,3]+integrate(function(k) otmFun(model.list.new,k=k,r=r,q=q,dT=t) * u * (u-1) * exp((u-1)*k) * exp(r*t),lower=L,upper=logF,subdivisions = 100L)$value})
  }

  if (doPlot>0) {
    t <- unique(regr.mat.copy$dT)[doPlot]
    pan.sub <- subset(regr.mat.copy,dT==t)
    k.pred <- seq(2*quantile(pan.sub$k,1e-4),2*quantile(pan.sub$k,1-1e-4),length.out=1000)
    otm.pred <- otmFun(model.list.new,k=k.pred,r=pan.sub$r[1],q=pan.sub$q[1],dT=t)
    layout(matrix(c(1,2),1,2))
    plot(k.pred,log(otm.pred),type="l")
    points(option.panels[[doPlot]]$k,log(option.panels[[doPlot]]$relMid),col="red") 
    plot(k.pred,otm.pred,type="l")
    points(option.panels[[doPlot]]$k,option.panels[[doPlot]]$relMid,col="red") 
    print(max(u.t.out[,3]))
  }
  return(u.t.out)
}