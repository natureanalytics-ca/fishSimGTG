#' Numerical simulations of fish population dynamics
#'
#' Conducts projections, MSE, and equilibrium conditions for a specified fish life history.
#'
#' @docType package
#' @name fishSimGTG
"_PACKAGE"


#' Life history example
#'
#' A LifeHistory object, prepopulated with an example fish life history.
#'
#'This life history is based on the Beverton-Holt life history invariants of M/K=1.5 and Lm/Linf=0.66
#' @format An S4 object of class LifeHistory
#' @source Simulated data
"LifeHistoryExample"

#' LBSPR sim example
#'
#' The data object LifeHistoryExample is used in lsbsprSimWrapper to produce YPR & SPR arrays.
#'
#' @format An S4 object of class LBSPRarray
#' @source Simulated data
"lbsprSimExample"

#' GTG sim example
#'
#' The data object LifeHistoryExample is used in gtgYPRWrapper to produce YPR & SPR arrays.
#'
#' @format An S4 object of class LBSPRarray
#' @source Simulated data
"gtgSimExample"
