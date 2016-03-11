# This script tests whether the Fourier-Transform  - Inverse Fourier Transform pairs hold in the data
library('affineOption',lib.loc=test.lib.dir)

# look at monthly horizon
for (dT in c(1/12,1)) {
print(dT)

params.BS <- namesToStruct(bsParamVec())

# create Heston model
params.BS$svFast$kpp <- 4
params.BS$svFast$lmb <- .2
params.BS$svFast$eta <- .04
params.BS$svFast$rho <- -.8

params.BS$jmp$lvec[1] <- 5
params.BS$jmp$muY <- -.15
params.BS$jmp$sigmaY <- .05

kVec <- seq(-16,8,by=.05)
strikeMat <- array(kronecker(sqrt(stateVec[,1])*sqrt(dT),kVec),c(1,length(kVec),nrow(stateVec)))
stateVec <- cbind(seq(.02,.06,.02),0)
colnames(stateVec) <- c("v1","v2")

odeSolFun <- function(x) {
  return(twoFactorJumpODEsSolve(cbind(1i*x,0,0),params.BS$svFast,params.BS$svSlow,params.BS$jmp,mkt=data.frame(t=dT,p=1,r=0,q=0),rtol=1e-12,atol=1e-12))
}
prices.jump<- gaussLaguerrePricer(odeSolFun,strikeMat=strikeMat,mkt=data.frame(t=dT,p=1,r=0,q=0),stateVec=stateVec,N=128)

uVec <- seq(-4,4)
for (uu in uVec) {
  ode.sol <- twoFactorJumpODEsSolve(u=matrix(c(uu,0,0),1,3),params.BS$svFast,params.BS$svSlow,params.BS$jmp,mkt=data.frame(t=dT,p=1,r=0,q=0))
  for (ii in 1:nrow(stateVec)) {
    # retain only prices that are reliable
    reliable.prices <- prices.jump$otm[1,,ii] > 1e-7
    igft.data <- impliedGFT(cbind(exp(strikeMat[1,reliable.prices,ii]),prices.jump$otm[1,reliable.prices,ii]),mkt=data.frame(t=dT,p=1,r=0,q=0),uu)
    igft.ode <- exp(ode.sol[1,1,"a"]+ode.sol[1,1,"b1"]*stateVec[ii,1]+ode.sol[1,1,"b2"]*stateVec[ii,2])
    expect_that(as.numeric(igft.data$cfI+igft.data$cfR),equals(as.numeric(igft.ode),tolerance=2e-4))
    #print(c(uu,ii,as.numeric(igft.data$cfI+igft.data$cfR),as.numeric(igft.ode)))
  }
}
}