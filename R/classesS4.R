

#----------------------
#Life history object
#----------------------

#Roxygen header
#'Life history object
#'
#'An S4 object that holds a description of a life history.
#'
#'This S4 object is used as input to a variety of fuctions, including various population dynamics simulation functions.
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
           ID = "character")
)

#lbspr Object
setClass("LBSPRarray",
  representation(
    LifeHistory = "LifeHistory",
    sim = "list"
  )
)

