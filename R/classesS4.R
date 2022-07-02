

#----------------------
#Life history object
#----------------------

#Roxygen header
#'Life history object
#'
#'An S4 object that holds a description of a life history.
#'
#'This S4 object is used as to input a life history.
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
#' @param L95delta Length increment between L50 and length at 95% maturity. Must be a value larger than 0.
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
#' @param H95delta Length increment between H50 and length 95% male. Must be a value larger than 0.
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
           L95delta = "numeric",
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
           H95delta = "numeric",
           author = "character",
           authAffiliation = "character",
           longDescription = "character",
           appBuild = "data.frame",
           ID = "character"
         ),
         prototype(
           title = "Example fish",
           speciesName = "Example fish",
           shortDescription = "Simulated life history of a fish based on B-H invariants",
           L_type = "TL",
           L_units = "cm",
           Walpha_units = "g",
           Linf = 100,
           K =  0.2,
           t0 = 0,
           L50 = 66,
           L95delta = 1,
           M =  0.3,
           MK = 1.5,
           LW_A = 0.01,
           LW_B = 3,
           Steep = 0.99,
           R0 = 1000,
           recSD = 0.6,
           recRho = 0,
           isHermaph = FALSE
         )
)


#----------------------
#Stock object
#----------------------

#Roxygen header
#'Fishery object
#'
#'An S4 object that holds a description of a fish stock, including selectivity and discard information.
#'
#'Options for vulnerability and retention retention functions along with guidance on parameter specification is found in Sel documentation
#' @param title A title for the object, useful for displaying the contents of the object
#' @param vulType String. Vulnerability function for historical time period, see selWrapper for options
#' @param vulParams Numeric value or vector for vulnerability params for historical time period. See selWrapper for options
#' @param retType String. Retention function for historical time period. See selWrapper for options
#' @param retParams Numeric value or vector for retention params for historical time period. See selWrapper for options
#' @param retMax Numeric value that defines the peak of the historical retention curve. A value between 0 and 1.
#' @param Dmort Discard mortality rate (not instantaneous rate, rather it is the fraction of discards killed e.g. 0.25 is 25% killed). A value between 0 and 1.
#' @importFrom methods new
#'
setClass("Fishery",
         representation(
           title = "character",
           vulType = "character",
           vulParams = "numeric",
           retType = "character",
           retParams = "numeric",
           retMax = "numeric",
           Dmort = "numeric"
         ),
         prototype (
           title = "Fishery corresponding to example fish",
           vulType = "logistic",
           vulParams = c(50,75),
           retType = "full",
           retMax = 1,
           Dmort = 0
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
#'Inputs for number of gtg, time step, and areas
#' @param title A title for the object, useful for displaying the contents of the object
#' @param gtg Number of growth-type groups
#' @param areas Number of areas in the model, must be greater than 1.
#' @param recArea A vector of length areas. Fraction of recruitment to each area with values summing to 1.
#' @param move A matrix of migration rates of dimensions areas x areas
#' @param iterations Number of iterations to run
#' @param historicalYears Number of years to simulate historical dynamics
#' @param historicalBio Number greater than 0 and less than 1. Model assumes we are dealing with an already exploited fish population
#' @param historicalBioType String. The type of historical biomass state, options are: 'relB' or 'SPR'.
#' @param historicalEffort A matrix of nrows = historicalYears and ncols = areas that contains value multipiers of initial equilibrium fishing effort
#' @importFrom methods new

setClass("TimeArea",
         representation(
           title = "character",
           gtg = "numeric",
           areas = "numeric",
           recArea = "numeric",
           move = "matrix",
           iterations = "numeric",
           historicalYears = "numeric",
           historicalBio = "numeric",
           historicalBioType = "character",
           historicalEffort = "matrix"
         ),
         prototype(
           title = "TimeArea corresponding to example fish",
           gtg = 13,
           areas = 2,
           recArea = c(0.99, 0.01),
           move = matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE),
           iterations = 5,
           historicalYears = 50,
           historicalBio = 0.5,
           historicalBioType = "relB",
           historicalEffort = matrix(1:1, nrow=50, ncol=2)
         )
)

#----------------------
#Strategy object
#----------------------
#Roxygen header
#'Strategy object
#'
#'An S4 object that holds descriptions of projections to be made, including a harvest strategy
#'
#'Details of projection to be made
#' @param title A title for the object, useful for describing the strategy
#' @param projectionYears Number of forward projection years to simulate
#' @param projectionName String. The name of projection method to apply. This is the name of the projection function
#' @param projectionParams List. List structure follows specification of the projection function specified in projectionName
#' @importFrom methods new

setClass("Strategy",
         representation(
            title = "character",
            projectionYears = "numeric",
            projectionName = "character",
            projectionParams = "list"
         )
)

#----------------------
#Stochastic object
#----------------------
#Roxygen header
#'Stochastic object
#'
#'An S4 object that holds parameters for stochastic components of the population dynamics.
#'
#'Details of stochastic components of the population dynamics. This list creates additional inputs as well as overrides for parameters specified elsewhere, allowing corresponding model components to become stochastic. Exception is recruitment variation, which is entered in the LifeHistory object
#' @param title A title for the object, useful for describing the scenario under exploration
#' @param historicalBio A vector of length 2 that contains a min and a max for historical equilibrium biomass. Replaces TimeArea@historicalBio. Continues to rely on TimeArea@historicalBioType. Range sampled at each iteration using a uniform distribution.
#' @param historicalCPUE A vector of length 2 that contains a min and a max for historical equilibrium fishery CPUE. Used to calculate a scaling coefficient, q, which allows results to be scaled to observed levels of CPUE. Range sampled at each iteration using a uniform distribution.
#' @param historicalCPUEType String. The type of fishery dependent CPUE, options are: 'vulB' or 'vulN'. Choice is used to calculating scaling coefficient relative to vulnerable biomass or vulnerable abundance.
#' @importFrom methods new

setClass("Stochastic",
         representation(
           title = "character",
           historicalBio = "numeric",
           historicalCPUE = 'numeric',
           historicalCPUEType = "character",
           Linf = "numeric",
           K = "numeric",
           L50 = "numeric",
           L95delta = "numeric",
           M = "numeric",
           Steep = "numeric",
           recSD = "numeric",
           recRho = "numeric",
           H50 = "numeric",
           H95delta = "numeric",
           histFisheryVul = "matrix",
           proFisheryVul = "matrix",
           sameFisheryVul = "logical",
           histFisheryRet = "matrix",
           proFisheryRet = "matrix",
           sameFisheryRet = "logical",
           histFisheryDmort = "matrix",
           proFisheryDmort = "matrix",
           sameFisheryDmort = "logical"
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
    lh = "list",
    sim = "list"
  )
)

#----------------------
#LBSRR sim ypr object
#----------------------

#Roxygen header
#'LBSPR YPR object
#'
#'An S4 object that holds the output of YPR analysis in a standardized format.
#'#' @importFrom methods new
setClass("LBSPRarray",
         representation(
           LifeHistory = "LifeHistory",
           sim = "list"
         )
)

