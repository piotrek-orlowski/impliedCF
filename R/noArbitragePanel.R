#' @title Check no-arbitrate conditions in an option panel
#' @description Given an option panel, the value of the underlying, the interest
#'   rate for the appropriate maturity and the dividend yield, check 
#'   no-arbitrage conditions for static payoffs, i.e. that there can be no
#'   non-positive value portfolio with a strictly positive payoff. This follows
#'   \code{http://www.fabiomercurio.it/NoArbitrage.pdf}
#' @param params.Q A list containing the Q parameter structures
#' @param stateMat An Sx2 matrix (with named columns v1,v2) that contains the 
#'   states at which the VIX is to be calculated
#' @param dT A positive scalar, which is the interval over which the VIX should 
#'   be calculated
#' @param panel data.frame of option prices with fields r (interest rate), logF (log-forward rate), dT (maturity), k (log-strike), relMid (relative price)
#' @export
#' @return is.no.arbitrage logical, 1 if there is no arbitrage

noArbitrageOptionCheck <- function(pn,r,q,s0,call=TRUE){
  
  ### create A matrix of strike transformations
  A <- matrix(0,nrow(pn)+2,nrow(pn)+2)
  A[,1] <- c(rep(1,nrow(pn)+1),0)
  A[nrow(A),] <- c(0,rep(1,nrow(pn)+1))
  A[-c(1,nrow(A)),2] <- exp(pn$k)
  for(nn in 3:ncol(A)){
    A[-c(1,nrow(A)),nn] <- exp(pn$k) - exp(pn$k[nn-2])
  }
  ### fill upper zeros
  diag.zeros <- upper.tri(A, diag = FALSE)
  A[which(diag.zeros)] <- 0
  
  ### create asset value vector
  T <- median(unique(pn$dT))
  v.bar <- c(exp(-r*T),s0*(exp(-q*T)),pn$relMid)
  
  ### Check
  
  is.no.arbitrage <- all(v.bar %*% solve(A) > 0)
  wrong.opts <- which(v.bar %*% solve(A) <= 0 + 1e-12)
  if(length(which(wrong.opts %in% c(1,2)))){
    wrong.opts <- wrong.opts[-which(wrong.opts %in% c(1,2))]
  }
  wrong.opts <- unique(wrong.opts) - 1
  if(length(wrong.opts) == 0){
    is.no.arbitrage <- 1
  }
  return(list(no.arb = is.no.arbitrage, wrong.opts = wrong.opts))
}

#' @export

removeArbitrageInstruments <- function(panel){
  
  r <- unique(pn$r)
  s0 <- 1
  logF <- unique(pn$logF)
  TTM <- median(unique(pn$dT))
  q <- r - logF/TTM
  
  pn.call <- pn[which(exp(pn$k) >= s0*exp((r-q)*TTM)),]
  pn.put <- pn[which(exp(pn$k) < s0*exp((r-q)*TTM)),]
 
  # change put prices to call prices via the put-call parity
  pn.put$relMid <- pn.put$relMid + s0*exp(q*TTM) - exp(pn.put$k)*exp(-r*TTM)
#   pn.call <- rbind(pn.put,pn.call)
  
  # remove bad call options
  arb.check.call <- noArbitrageOptionCheck(pn.call,r,q,s0,TRUE)
  wrong.opts.list <- list()
  while(!arb.check.call$no.arb){
    
    pn.call <- pn.call[-arb.check.call$wrong.opts,]
    arb.check.call <- noArbitrageOptionCheck(pn.call,r,q,s0,TRUE)
  }
  
  # remove bad put options
  
  arb.check.put <- noArbitrageOptionCheck(pn.put,r,q,s0,TRUE)
  iter <- 1
  while(!arb.check.put$no.arb & !identical(matrix(1),matrix(arb.check.put$wrong.opts))){
    pn.put <- pn.put[-arb.check.put$wrong.opts,]
    arb.check.put <- noArbitrageOptionCheck(pn.put,r,q,s0,TRUE)
    iter <- iter + 1
  }
  # back via the pc-parity
  pn.put$relMid <- pn.put$relMid + exp(pn.put$k)*exp(-r*TTM) - s0*exp(q*TTM)

  pn.call <- rbind(pn.put,pn.call)

  # take put prices back to be puts
  return(pn.call)
}

removeArbitrageJointly <- function(pn, r, q, s0){
  
  TTM <- median(unique(pn$dT))
  pn.call <- pn[which(exp(pn$k) >= s0*exp((r-q)*TTM)),]
  pn.put <- pn[which(exp(pn$k) < s0*exp((r-q)*TTM)),]
  
  # change put prices to call prices via the put-call parity
  pn.put$relMid <- pn.put$relMid + s0*exp(q*TTM) - exp(pn.put$k)*exp(-r*TTM)
  pn.call <- rbind(pn.put,pn.call)
  
  # remove bad options
  arb.check.call <- noArbitrageOptionCheck(pn.call,r,q,s0,TRUE)
  wrong.opts.list <- list()
  while(!arb.check.call$no.arb){
    
    pn.call <- pn.call[-arb.check.call$wrong.opts,]
    arb.check.call <- noArbitrageOptionCheck(pn.call,r,q,s0,TRUE)
  }
  # back via the pc-parity
  pn.call$relMid[which(exp(pn.call$k) < unique(pn.call$logF))] <- pn.call$relMid[which(exp(pn.call$k) < unique(pn.call$logF))] + exp(pn.call$k[which(exp(pn.call$k) < unique(pn.call$logF))])*exp(-r*TTM) - s0*exp(q*TTM)
  
  return(pn.call)
}