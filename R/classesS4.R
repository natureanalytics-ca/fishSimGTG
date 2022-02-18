

#----------------------
#Life history object
#----------------------

#Roxygen header
#'Life history object
#'
#'An S4 object that holds a description of a life history.
#'
#'This S4 object is used as input to a variety of functions, including various population dynamics simulation functions.
#' @param title A title for the object, useful for displaying the contents of the object
#' @param speciesName Scientific name of the species
#' @param shortDescription A brief description of the object. This could be the common name, stock, geographic location of the stock, etc.
#' @param L_type The method of measuring length. e.g. TL for total length. Must be consistent for all length params e.g., Linf, L50, L95
#' @param L_units Units of measure for the object. cm is expected. Must be consistent for all length params e.g., Linf, L50, L95
#' @param W_units Units of weight for the object. Must be consistent with LW_A and LW_B params
#' @param Linf von Bertalanffy Loo parameter
#' @param K von Bertalanffy K parameter per year
#' @param t0 von Bertalanffy t0 parameter
#' @param L50 Length at 50% maturity.
#' @param L95 Length at 95% maturity. Must be a value larger than L50
#' @param M Natural mortality rate per year
#' @param MK Ratio of M to K
#' @param LW_A Parameter for length-weight relationship W=aL^b
#' @param LW_B Parameter for length-weight relationship W=aL^b
#' @param Tmax Maximum observed age
#' @param Steep Steepness of the Beverton-Holt stock recruit relationship
#' @param R0 Unfished number of recruits
#' @param recSD Inter annual recruitment variation
#' @param recRho Inter-annual correlation in recruitment variation
#' @param isHermaph Logical whether species is a protogynous hermaphrodite (TRUE). FALSE is gonochoristic species.
#' @param H50 Length at 50% male.
#' @param H95 Length at 95% male. Must be a value larger than H50
#' @param author Author of the life history - do not include your email address as these life histories may be posted online
#' @param authAffiliation A way to identify the author without relying on an email address
#' @param longDescription Document the rational for choices made in creating the life history. Other users will rely on this information.
#' @param appBuild A data frame that holds details when the life history is built using the Shiny app
#' @param ID A user ID. Not required.
#' @importFrom methods new


setClass("LifeHistory",
         representation(
           title = "character",
           speciesName = "character",
           shortDescription = "character",
           L_type = "character",
           L_units = "character",
           Walpha_units = "character",
           Linf = "numeric",
           K =  "numeric",
           t0 = "numeric",
           L50 = "numeric",
           L95 = "numeric",
           M =  "numeric",
           MK = "numeric",
           LW_A = "numeric",
           LW_B = "numeric",
           Tmax = "numeric",
           Steep = "numeric",
           R0 = "numeric",
           recSD = "numeric",
           recRho = "numeric",
           isHermaph = "logical",
           H50 = "numeric",
           H95 = "numeric",
           author = "character",
           authAffiliation = "character",
           longDescription = "character",
           appBuild = "data.frame",
           ID = "character")
)


#----------------------
#Stock object
#----------------------

#Roxygen header
#'Stock object
#'
#'An S4 object that holds a description of a fish stock, including selectivity and discard information.
#'
#'This S4 object is used as input to a variety of functions, including various population dynamics simulation functions.
#'Options for vulnerability and retention retention functions along with guidance on parameter specification is found in Sel documentation
#' @param title A title for the object, useful for displaying the contents of the object
#' @param historicalVulType String. Vulnerability function for historical time period, see selWrapper for options
#' @param historicalVulParams Numeric value or vector for vulnerability params for historical time period. See selWrapper for options
#' @param historicalRetType String. Retention function for historical time period. See selWrapper for options
#' @param historicalRetParams Numeric value or vector for retention params for historical time period. See selWrapper for options
#' @param historicalRetMax Numeric value that defines the peak of the historical retention curve. A value between 0 and 1.
#' @param historicalDmort Historical tiem period discard mortality rate (not instantaneous rate, rather it is the fraction of discards killed e.g. 0.25 is 25% killed). A value between 0 and 1.
#' @param projectionVulType String. Vulnerability function for projection time period, see selWrapper for options
#' @param projectionVulParams Numeric value or vector for vulnerability params for projection time period. See selWrapper for options
#' @param projectionRetType String. Retention function for projection time period. See selWrapper for options
#' @param projectionRetParams Numeric value or vector for retention params for projection time period. See selWrapper for options
#' @param projectionRetMax Numeric value that defines the peak of the projection retention curve. A value between 0 and 1.
#' @param projectionDmort Projection time period discard mortality rate (not instantaneous rate, rather it is the fraction of discards killed e.g. 0.25 is 25% killed). A value between 0 and 1.
#' @importFrom methods new
#'
setClass("Fishery",
         representation(
           title = "character",
           historicalVulType = "character",
           historicalVulParams = "numeric",
           historicalRetType = "character",
           historicalRetParams = "numeric",
           historicalRetMax = "numeric",
           historicalDmort = "numeric",
           projectionVulType = "character",
           projectionVulParams = "numeric",
           projectionRetType = "character",
           projectionRetParams = "numeric",
           projectionRetMax = "numeric",
           projectionDmort = "numeric"
         )
)

#----------------------
#Time-area object
#----------------------
#Roxygen header
#'Time-area object
#'
#'An S4 object that holds descriptions of time step, gtg, and area params
#'
#'This S4 object is used as input to a variety of functions, including various population dynamics simulation functions.
#'Options for vulnerability and retention retention functions along with guidance on parameter specification is found in Sel documentation
#' @param title A title for the object, useful for displaying the contents of the object
#' @param gtg Number of growth-type groups
#' @param stepsPerYear Number of time steps per year (e.g., 1 for annual time step, 12 for monthly time step)
#' @importFrom methods new

setClass("TimeArea",
         representation(
           title = "character",
           gtg = "numeric",
           stepsPerYear = "numeric"
         )
)

#----------------------
#YPR object
#----------------------

#Roxygen header
#'YPR object
#'
#'An S4 object that holds the output of YPR analysis in a standardized format.
#'#' @importFrom methods new
setClass("YPRarray",
  representation(
    LifeHistory = "LifeHistory",
    sim = "list"
  )
)


