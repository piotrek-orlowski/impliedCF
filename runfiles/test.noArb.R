#' Test the no-arbitrage filter based on quantile regression
rm(list=ls())
load("D://optionPanels_2013.RData")

# set day
i <- 14630

plot(optionPanelsOTM[[i]]$k,optionPanelsOTM[[i]]$relMid)

option.noarb <- noArbPrices(optionPanelsOTM[[i]])

lines(option.noarb$k,option.noarb$relMid)

plot(option.noarb$k,option.noarb$relMid/optionPanelsOTM[[i]]$relMid-1,type="l")
