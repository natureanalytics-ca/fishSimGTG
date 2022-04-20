


#----------------------------------------------------------------------
#Fixed time period - for fixed harvest time period (not for MSE time period)
#----------------------------------------------------------------------

#Roxygen header
#'Historical fishing pressure
#'
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @export

fixedStrategy<-function(phase, dataObject){

  #Unpack dataObject
  j <- areas <- k <- TimeAreaObj <- is <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #Booking keeping for year for items in TimeAreaObj
  yr <- j - 1

  if(phase==3){
    #Create a temp data frame of fishing mortalities by area
    Flocal<-data.frame()
    for (m in 1:areas) Flocal<-rbind(Flocal, c(j, k, m, TimeAreaObj@historicalEffort[yr,m]*is$Feq))
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}


#-------------------------------------------------------------------
#Projection modeling - no harvest control rule, simple projections
#-------------------------------------------------------------------

#Roxygen header
#' Static projections of combinations of effort changes, bag limit, and/or spatial closures. Also used in combination with ProFisheryObj (e.g. size limit change) to project temporal dynamics of size limits.
#'
#' The Strategy object should be specified as follows. (1) Set Strategy@projectionYears to the number of forward projection years you wish to simulate. (2) Strategy@projectionName = "projectionStrategy".
#' (3)  Strategy@projectionParams should be a list with two items. First items is a vector of length areas containing bag limit. For no bag limit use -99.
#' The bag limit should be thought of as take per unit time (e.g. day) and basically acts like a CPUE threshold.
#' The effect of bag limit is calculated against the historical CPUE (e.g. in same units of take per unit time) in the Stochastic object historicalCPUE. Make sure that the bag limit and historicalCPUE are consistent with historicalCPUEType (e.g., biomass or abundance (numbers)) and this parameter is used in determing the effect of fish biomass or abundance on CPUE.
#' The second item in the list is a matrix of nrows = projectionYears and ncols = areas that contains value multipiers of initial equilibrium fishing effort. This allows projection of effort reduction and of marine reserves via setting effort to 0.
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @importFrom stats dpois ppois
#' @export

projectionStrategy<-function(phase, dataObject){

  #Unpack dataObject
  j <- TimeAreaObj <- areas <- StrategyObj <- is <- k <- StochasticObj <- lh <- N <- selHist <- Cdev <- selGroup <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #Book keeping year for items in StrategyObj
  yr <- j - TimeAreaObj@historicalYears - 1

  #Book keeping year for terminal year of historical period, if available
  yrHist <- TimeAreaObj@historicalYears + 1

  if(phase==3){

    Flocal<-data.frame()
    for (m in 1:areas){

      bag <- StrategyObj@projectionParams[[1]][m]

      if(bag == -99){

        #Apply Flocal
        Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*Ftotal[yrHist,k,m]
        Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

      } else {

        if(StochasticObj@historicalCPUEType == "vulN") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist$keep[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup$keep[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*Ftotal[yrHist,k,m]*Fmult
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }

        if(StochasticObj@historicalCPUEType == "vulB") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist$keep[[x]]*lh$W[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup$keep[[x]]*lh$W[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*Ftotal[yrHist,k,m]*Fmult
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }
      }
    }
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}
