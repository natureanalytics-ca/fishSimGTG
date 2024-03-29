% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YPR.R
\name{gtgYPRWrapper_Fonly}
\alias{gtgYPRWrapper_Fonly}
\title{YPR simulation using growth-type group model for sequence of F_M. Based on pre-specified LHWrapper and selWrapper objects.}
\usage{
gtgYPRWrapper_Fonly(lh, sel, F_Step = 0.02)
}
\arguments{
\item{lh}{An object produced by LHWrapper.}

\item{sel}{An object produced by selWrapper}

\item{F_Step}{F step size for sequence of F_M}
}
\description{
Produces YPR, Yield, and SPR arrays across combinations of F_M
}
\details{
Utilizes the growth-type group model to calculate YPR, Yield, and SPR across sequency of F_M.
Produces quantities of interest like MSY, etc.
Required parameters of the LifeHistory object are: Linf, L50, L95delta, M, K, Steep, LW_A, LW_B
}
