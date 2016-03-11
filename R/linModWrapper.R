#' Wrapper to \code{\link{solve.QP}} for estimating constrained linear models.
#' 
#' @aliases constr.lm.estim
#' @aliases pred.impliedCF.lm
#' @title Estimation of constrained linear model, prediction.
#' @description This function takes a formula, a vector of term formulae and estimates a constrained linear model, giving it class \code{impliedCF.lm} so that the predict routine can be used as well.
#' @param terms Formulae for left and right hand side terms (individual components of formula). The explained variables is assumed to be in \code{terms[1]}
#' @param data Data frame with right hand side variables
#' @param R Matrix of linear restrictions of the form Rb >= r
#' @param r Right hand side of linear restriction Rb >= r
#' @param weights For some weighted least squares
#' @param meq check \code{\link{quadprog::solve.QP}}
#' @rdname linModWrapper
#' @export
#' 
constr.lm.estim <- function(terms, data, R, r, weights = rep(1,nrow(data)), meq = 0){
  
  stopifnot(!is.null(terms))
  stopifnot(!is.null(data))
  
  # Evaluate terms from data
  Xmat <- vapply(X = terms[2:length(terms)], FUN = function(ter){
    return(eval(parse(text=ter), envir = data))  
  }, FUN.VALUE = rep(0,nrow(data)))
  Xmat.orig <- Xmat
  
  # Obtain explained variables
  Yvec <- matrix(eval(parse(text = terms[1]), envir = data), nrow(data), 1)
  
  # Multiply by weights
  Yvec <- weights * Yvec
  Xmat <- apply(X = Xmat, MARGIN = 2, FUN = '*', weights)
  
  # Create matrices that solve.QP understands
  Dmat <- t(Xmat) %*% Xmat
  if(min(eigen(Dmat)$values) < 0){
    Dmat <- makePsd(Dmat)
  }
  while(min(eigen(Dmat)$values) < 0){
    Dmat <- makePsd(Dmat)
  }
  dvec <- t(Yvec) %*% Xmat
  
  # Solve
  lm <- solve.QP(Dmat, dvec, t(R), r, factorized = FALSE, meq)
  
  # Make solution structure of a certain class
  sol <- list()
  sol$coefficients <- lm$solution
  sol$rss <- lm$value
  sol$unc <- lm$unconstrained.solution
  sol$act <- lm$iact
  sol$fitted.values <- Xmat.orig %*% matrix(sol$coefficients, nrow = length(sol$coefficients), ncol = 1)
  sol$terms <- terms
  
  class(sol) <- c("impliedCF.lm")
  
  return(sol)
}

#' @rdname linModWrapper
#' @param model model of class \code{impliedCF.lm}
#' @param newdata conforming data frame
#' @export
#' 
predict.impliedCF.lm <- function(model, newdata){
  
  terms <- model$terms
  # Evaluate terms from data
  start.at <- 1
  if(length(terms) != length(model$coefficients)){
    start.at <- 2
  }
  Xmat <- vapply(X = terms[start.at:length(terms)], FUN = function(ter){
    return(eval(parse(text=ter), envir = newdata))  
  }, FUN.VALUE = rep(0,nrow(newdata)))
  
  beta.hat <- matrix(model$coefficients, nrow = length(model$coefficients), ncol = 1)
  
  # Calculate prediction
  pred <- Xmat %*% beta.hat
  
  return(pred)
}