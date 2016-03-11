#' Calculate the theoretical value of the VIX
#' 
#' This function calculates the theoretical VIX value, by differentiating the ODE solution around 0.
#' @param params.Q A list containing the Q parameter structures
#' @param stateMat An Sx2 matrix (with named columns v1,v2) that contains the states at which the VIX is to be calculated
#' @param dT A positive scalar, which is the interval over which the VIX should be calculated
#' @export
#' @return A vecotr of length S, with the theoretical VIX values
VIXode <- function(params.Q,stateMat,dT) {
  eps <- 1e-6
  j <- 1
  uVec <- eps * (1:3)
  ode.deriv <- twoFactorJumpODEsSolve(cbind(uVec,0,0),params.Q$svFast,params.Q$svSlow,params.Q$jmp,mkt=list(p=1,r=0,q=0,t=dT))
  # calculate the first derivative, which is the first cumulant
  var.coeffs <- rep(NA,3)
  for (ii in 1:3) {
    var.coeffs[ii] <--2* diff(ode.deriv[,1,ii],differences=j)[1]*eps^-j /dT
  }
  vix <- as.numeric(var.coeffs[3]+var.coeffs[1]*stateMat[,1]+var.coeffs[2]*stateMat[,2])
  return(vix)
}