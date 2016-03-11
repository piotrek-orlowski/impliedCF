#' @title Implied divergences from option DB file
#' @description This function takes a data.frame of option data and at each date builds the spline IV surface, based on that estimates prices of divergence, skewness and quarticity swaps
#' @param options data.frame with fields date, k (log(relative strike with respect to underlying)), relMid (relative mid-price with respect to underlying), r (interest rate), q (implied dividend yield), k_stand (k divided by atmIV and sqrt(expiration))
#' @param pow.vec vector of divergence powers
#' @param t.vec vector of maturities
#' @param types vector of characters, length up to 3, can contain "div", "skew" and "quart"
#' @param filtering.expression quoted expression to be passed to \code{\link[dplyr]{filter_}}, to select all options take something that evaluates to TRUE, like \code{quote(TRUE)}, to select options by strike, for example \code{quote(k_stand <= 3 & k_stand >= -6)}
#' @param cl cluster for parallel calculations
#' @param U upper integration limit for pricing in terms of log-strike. U > 0
#' @param L lower integration limit for pricing in terms of log-strike. L < 0
#' @param gam.bs bs argument to \code{\link[mgcv]{s}}
#' @param gam.m m argument to \code{\link[mgcv]{s}}
#' @export

impliedDivergencesWrapper <- function(options, pow.vec, t.vec, types, filtering.expression, cl = cl, U = NULL, L = NULL, gam.bs = "ds", gam.m = c(1,0.5)){
  
  ut.mat <- expand.grid(u=pow.vec, t=t.vec,type = types)
  
  options <- options %>% filter_(filtering.expression)
  
  option.panels <- vector(mode = "list", length = length(unique(options$date)))
  for(ind in 1:length(option.panels)){
    un.date <- unique(options$date)[ind]
    loc.opts <- options %>% filter(date == un.date)
    expiries <- unique(loc.opts$expiration)
    loc.list <- vector(mode = 'list', length = length(expiries))
    loc.mkt <- data.frame(p = numeric(), q = numeric(), r = numeric(), t = numeric())
    for(ep in 1:length(expiries)){
      loc.loc.opts <- loc.opts %>% filter(expiration == expiries[ep])
      loc.list[[ep]] <- loc.loc.opts %>% ungroup %>% select(k,relMid) %>% as.matrix
      loc.mkt <- rbind(loc.mkt, data.frame(p=1,r=unique(loc.loc.opts %>% .$r), t= unique(loc.loc.opts %>% .$expiration), q = unique(loc.loc.opts %>% .$q)))
    }
    option.panels[[ind]] <- list(opt.pn = loc.list, mkt = loc.mkt, day = un.date)
  }
  
  if(!is.null(cl)){
    clusterExport(cl,c("pow.vec","t.vec","types","U","L"))
    clusterEvalQ(cl, ut.mat <- expand.grid(u=pow.vec, t=t.vec,type = types))
    
    div.pr.db <- parLapplyLB(cl = cl, X = option.panels, fun = function(opts){
      res <- tryCatch(expr = tpsImpliedDivergence(option.panels = opts$opt.pn, mkt.frame = opts$mkt, u.t.mat = ut.mat, verbose = F, time.IV = T, L = L, U = U),
                      error = function(e){
                        print(e)
                        return(data.frame(u=NA_real_, t = NA_real_, type = NA_character_, res = NA_real_, stringsAsFactors = F))
                      }
      )
      res <- cbind(day = opts$day, res)
      return(res)
    })
  } else {
    div.pr.db <- lapply(X = option.panels, FUN = function(opts){
      res <- tryCatch(expr = tpsImpliedDivergence(option.panels = opts$opt.pn, mkt.frame = opts$mkt, u.t.mat = ut.mat, verbose = F, time.IV = T, L = L, U = U,  gam.bs = gam.bs, gam.m = gam.m),
                      error = function(e){
                        print(e)
                        return(data.frame(u=NA_real_, t = NA_real_, type = NA_character_, res = NA_real_, stringsAsFactors = F))
                      }
      )
      res <- cbind(day = opts$day, res)
      return(res)
    })
  }
  
  div.pr.db <- rbind_all(div.pr.db)
  return(div.pr.db)
}