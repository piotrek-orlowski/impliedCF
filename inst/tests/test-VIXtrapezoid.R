# create Heston prices and use it to test the pricer & VIX calculator
library('affineOption',lib.loc=test.lib.dir)

dT <- 1/2

params.BS <- namesToStruct(bsParamVec())

# create Heston model
params.BS$svFast$kpp <- 4
params.BS$svFast$lmb <- .2
params.BS$svFast$eta <- .04
params.BS$svFast$rho <- -.8

# create cf function for pricing
odeSolFun <- function(x) {
  return(twoFactorJumpODEsSolve(cbind(1i*x,0,0),params.BS$svFast,params.BS$svSlow,params.BS$jmp,mkt=data.frame(t=dT,p=1,r=0,q=0)))
}

# create strike matrix and initial states
kVec <- seq(-4,3,by=.1)
strikeMat <- array(kronecker(sqrt(stateVec[,1])*sqrt(dT),kVec),c(1,length(kVec),nrow(stateVec)))
stateVec <- cbind(seq(.02,.06,.02),0)
colnames(stateVec) <- c("v1","v2")

# calculate prices
prices.heston <- gaussLaguerrePricer(odeSolFun,strikeMat=strikeMat,mkt=data.frame(t=dT,p=1,r=0,q=0),stateVec=stateVec)

# now calculate VIX coefficient loadings (how they depend on initial volatility)
vix.ode <- VIXode(params.BS,stateVec,dT)

# check if theoretical and empirical vix values are aligned
for (nn in 1:nrow(stateVec)) {
  priceMat <- cbind(exp(strikeMat[,,nn]),prices.heston$otm[,,nn])
  vix.emp <- VIXtrapezoid(priceMat,mkt=list(t=dT,p=1,r=0,q=0))
  expect_that(vix.emp,equals(vix.ode[nn],tolerance=5e-3))
}

# now compare with closed form heston
for (nn in 1:nrow(stateVec)) {
  vix.heston <- (stateVec[nn,1]-params.BS$svFast$eta)/params.BS$svFast$kpp*(1-exp(-params.BS$svFast$kpp*dT))/dT + params.BS$svFast$eta
  expect_that(as.numeric(vix.heston),equals(vix.ode[nn],tolerance=5e-3))
}

############ now do the same exercise with jump
params.BS$jmp$lvec[1] <- 1
params.BS$jmp$muY <- -.15
params.BS$jmp$sigmaY <- .05

# create cf function for pricing
odeSolFun <- function(x) {
  return(twoFactorJumpODEsSolve(cbind(1i*x,0,0),params.BS$svFast,params.BS$svSlow,params.BS$jmp,mkt=data.frame(t=dT,p=1,r=0,q=0)))
}

# create strike matrix and initial states
kVec <- seq(-3.5,2.5,by=.1)
strikeMat <- array(kronecker(sqrt(stateVec[,1])*sqrt(dT),kVec),c(1,length(kVec),nrow(stateVec)))
stateVec <- cbind(seq(.02,.06,.02),0)
colnames(stateVec) <- c("v1","v2")

# calculate prices
prices.jump <- gaussLaguerrePricer(odeSolFun,strikeMat=strikeMat,mkt=data.frame(t=dT,p=1,r=0,q=0),stateVec=stateVec)

# now calculate VIX coefficient loadings (how they depend on initial volatility)
vix.ode <- VIXode(params.BS,stateVec,dT)

# check if theoretical and empirical vix values are aligned
for (nn in 1:nrow(stateVec)) {
  priceMat <- cbind(exp(strikeMat[,,nn]),prices.jump$otm[,,nn])
  vix.emp <- VIXtrapezoid(priceMat,mkt=list(t=dT,p=1,r=0,q=0))
  expect_that(vix.emp,equals(vix.ode[nn],tolerance=2e-2))
}