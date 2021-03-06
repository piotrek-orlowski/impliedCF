% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optionPanelInterpolate.R
\name{optionPanelInterpolate}
\alias{optionPanelInterpolate}
\title{Fit a smooth IV surface to all available options and return option prices for given maturities and strike ranges that depend on maturity and IV.}
\usage{
optionPanelInterpolate(option.panels, mkt.frame, out.mat, doPlot = FALSE,
  doFitPlot = 0, verbose = FALSE)
}
\arguments{
\item{option.panels}{A list, each element containing a dataframe with normalized mid OTM option prices and LOG strikes. Each panel should correspond to the SAME day.}

\item{out.mat}{Sx2 matrix, where first column specifies target maturities and the second specifies the number of options available at each maturity (for now they have to be equal)}

\item{doPlot}{Whether fitted / true values should be plotted}

\item{verbose}{Do we want diagnostic information?}

\item{mkt.list}{A data-frame, where each row corresponds to the market features associated with the n-th option panel (e.g. maturity, interest rate, etc...)}
}
\value{
out.panels list with fields \code{p.list} -- a list of option panels, \code{m.list} -- a list of corresponding market structures.
}
\description{
Fit a smooth IV surface to all available options and return option prices for given maturities and strike ranges that depend on maturity and IV.
}
