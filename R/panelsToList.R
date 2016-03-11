#' Convert joined panel into seperate lists, that can be fed to the igft interpolator
#' 
#' @param panel.names Vector of panel names that will be loaded and parsed
#' @export
#' @return Returns a list with components p.list, m.list and days.all. These correspond daily priceMat panel lists, daily market frames and the vector of available days.
panelsToList <- function(panel.names, noisy=TRUE) {
  # first load all the panels and create a list of priceMats, mkts and available maturities
  n.cols <- 0
  n.rows <- 0
  mkt.all.list <- priceMat.day.list <- priceMat.all.list <- list()
  
  curr.day <- -1
  curr.q <- curr.r <- curr.mat <- 0
  
  priceMat <- NULL
  mkt.frame <- NULL
  n.days <- 0
  
  days.all <- NULL
  
  for (pn in panel.names) {
    var.names <- load(pn)
    stopifnot("loc.opts" %in% var.names)
    for (oo in 1:(nrow(loc.opts)+1)) {
      # if a new day, save the list of panels and re-initialize the panel list
      if (loc.opts[oo,"day"] != curr.day || oo == nrow(loc.opts)+1) {
        if (!is.null(priceMat)) {
          # make sure last daily panel is saved
          n.panels <- length(priceMat.day.list)
          priceMat.day.list[[n.panels+1]] <- priceMat
          mkt.frame <- rbind(mkt.frame,c(curr.mat,curr.q,curr.r,1))
          
          # now save into the big list
          n.days <- n.days + 1
          priceMat.all.list[[n.days]] <- lapply(priceMat.day.list,function(x) {y <- data.frame(k=log(x[,1]),relMid=x[,2]);return(y)})
          priceMat.day.list <- list()
          colnames(mkt.frame) <- c("t","q","r","p")
          mkt.all.list[[n.days]] <- as.data.frame(mkt.frame)
        }
        priceMat <- NULL  
        mkt.frame <- NULL
        # stop if we are at the end
        if (oo == nrow(loc.opts)+1) {break}

        curr.day <- loc.opts[oo,"day"]
        curr.mat <- loc.opts[oo,"dT"]
        curr.r <- loc.opts[oo,"r"]
        curr.q <- loc.opts[oo,"q"]
        
        days.all <- c(days.all,curr.day)
      }
      # if new maturity, then save panel into daily list and re-initialize the daily priceMat
      if (loc.opts[oo,"dT"] != curr.mat & !is.null(priceMat)) {
        n.panels <- length(priceMat.day.list)
        priceMat.day.list[[n.panels+1]] <- priceMat
        priceMat <- NULL
        mkt.frame <- rbind(mkt.frame,c(curr.mat,curr.q,curr.r,1))
        
        curr.mat <- loc.opts[oo,"dT"]
        curr.r <- loc.opts[oo,"r"]
        curr.q <- loc.opts[oo,"q"]
      }
      if (noisy) {
        priceMat <- rbind(priceMat,c(exp(loc.opts[oo,"relStrike"]),loc.opts[oo,"otm.noise"]))
      } else {
        priceMat <- rbind(priceMat,c(exp(loc.opts[oo,"relStrike"]),loc.opts[oo,"otm"]))
      }
    }
  }
  return(list(p.list = priceMat.all.list, m.list = mkt.all.list, days.all = days.all))
}