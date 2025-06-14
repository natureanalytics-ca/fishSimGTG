% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classesS4.R
\docType{class}
\name{Stochastic-class}
\alias{Stochastic-class}
\title{Stochastic object}
\arguments{
\item{title}{A title for the object, useful for describing the scenario under exploration}

\item{historicalBio}{A vector of length 2 that contains a min and a max for historical equilibrium biomass. If entered, replaces TimeArea@historicalBio. Continues to rely on TimeArea@historicalBioType. Range sampled at each iteration using a uniform distribution.}

\item{Linf}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{K}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{L50}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{L95delta}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{M}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{Steep}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{recSD}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{recRho}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{H50}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{H95delta}{A vector of length 2 that contains a min and a max. If entered, values replaces value in life history object. Create unique value for each iteration by sampling from uniform distribution.}

\item{histFisheryVul}{A matrix n cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces HistFisheryObj@vulParams. Columns correspond to needed inputs of HistFisheryObj@vulType. Range sampled at each iteration using a uniform distribution.}

\item{proFisheryVul_list}{A list containing number of objects equal to areas. Each object is a matrix n cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces ProFisheryObj@vulParams. Columns correspond to needed inputs of ProFisheryObj@vulType. Range sampled at each iteration using a uniform distribution.}

\item{sameFisheryVul}{Logical. Indicates whether values generated for histFisheryVul should be applied so that historical and projection parameter values are identical. TRUE also overrides any input in proFisheryVul_list}

\item{histFisheryRet}{A matrix n cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces HistFisheryObj@retParams. Columns correspond to needed inputs of HistFisheryObj@retType. Range sampled at each iteration using a uniform distribution.}

\item{proFisheryRet_list}{A list containing number of objects equal to areas. Each object is a matrix n cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces ProFisheryObj@retParams. Columns correspond to needed inputs of ProFisheryObj@retType. Range sampled at each iteration using a uniform distribution.}

\item{sameFisheryRet}{Logical. Indicates whether values generated for histFisheryRet should be applied so that historical and projection parameter values are identical. TRUE also overrides any input in proFisheryRet_list}

\item{histFisheryDmort}{A matrix 1 cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces HistFisheryObj@Dmort. Range sampled at each iteration using a uniform distribution.}

\item{proFisheryDmort_list}{A list containing number of objects equal to areas. Each object is a matrix 1 cols and 2 rows, with rows 1 and 2 containing a min and a max for parameter corresponding to column n. If entered, replaces ProFisheryObj@Dmort. Range sampled at each iteration using a uniform distribution.}

\item{sameFisheryDmort}{Logical. Indicates whether values generated for histFisheryDmort should be applied so that historical and projection parameter values are identical. TRUE also overrides any input in proFisheryDmort_list}
}
\description{
An S4 object that holds parameters for stochastic components of the population dynamics.
}
\details{
Details of stochastic components of the population dynamics. This list creates additional inputs as well as overrides for parameters specified elsewhere, allowing corresponding model components to become stochastic. Exception is recruitment variation, which is entered in the LifeHistory object
}
