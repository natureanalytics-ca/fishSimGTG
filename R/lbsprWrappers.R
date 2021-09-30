

#----------------------
#LBSPR sim wrapper
#----------------------

#Roxygen header
#'LBSPR eumetric simulation
#'
#'Wrapper for LBSPRsim that produces YPR, Yield, and SPR arrays across combinations of F_M and Lc
#'
#'Utilizes the LBSPR library to calculate YPR, Yield, and SPR across factorial combination of F_M and Lc.
#'This produces an array of values for each output, which is useful for surface plot, finding MSY, etc.
#'Yield is presented in relative terms, with a maximum of 1. YPR is not presented in relative terms, but the user
#'may want to re-calcuate this quantity as relative YPR for presentation purposes.
#'Required parameters of the LifeHistory object are: Linf, L50, L95, MK, Steep, LW_A, LW_B
#' @param LifeHistory  A life history object.
#' @param binWidth LBSPR length bin width, default value is 1.
#' @param binMin LBSPR binMin, default is 0
#' @param progressName When used within a shiny app, this function can update a shinyWidgets::progress bar.
#' progressName is the id of the shinyWidgets::progress bar specified in a shiny app.
#' @param sessionName Required to update a shinyWidgets::progress bar
#' @import LBSPR
#' @importFrom shinyWidgets updateProgressBar
#' @importFrom methods new
#' @export
#' @examples
#' lh<-new("LifeHistory")
#' lh@Linf<-100
#' lh@L50<-66
#' lh@L95<-67
#' lh@MK<-1.5
#' lh@LW_A<-0.01
#' lh@LW_B<-3
#'
#' sim<-lbsprSimWrapper(lh)

lbsprSimWrapper<-function(LifeHistory, binWidth=1, binMin=0, progressName=NULL, sessionName=NULL){

  #-----------------------------
  #Create life history pars list
  #-----------------------------

  MyPars <- new("LB_pars")
  MyPars@Linf <- LifeHistory@Linf
  MyPars@L50 <- LifeHistory@L50
  MyPars@L95 <- LifeHistory@L95
  MyPars@MK <- LifeHistory@MK
  MyPars@BinWidth <- binWidth
  MyPars@Steepness<-LifeHistory@Steep
  MyPars@L_units <- LifeHistory@L_units
  MyPars@Walpha <- LifeHistory@LW_A
  MyPars@Walpha_units <-LifeHistory@Walpha_units
  MyPars@Wbeta <- LifeHistory@LW_B
  MyPars@FecB <- LifeHistory@LW_B
  MyPars@BinMin <- binMin
  #Setup place holder values for these parameters, we will change these later
  MyPars@SL50 <- LifeHistory@L50
  MyPars@SL95 <- LifeHistory@L50+1
  MyPars@FM<-1


  #------------------
  #Eumetric analysis
  #------------------

  Lmax<-(1 - 0.01^(1/MyPars@MK)) * MyPars@Linf
  Lc<-seq(floor(0.1*Lmax),  floor(Lmax), 1)
  F_M<-round(seq(0, 4, 0.2), 1)
  SPR_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))
  YPR_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))
  Yield_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))

  show_condition <- function(code) {
    tryCatch({
      x<-code
      list(SPR=x@SPR, YPR=x@YPR, Yield=x@Yield)
    },       error = function(c)NA,
    warning = function(c) NA,
    message = function(c) NA
    )
  }

  steps<-NROW(Lc)*NROW(F_M)
  counter <- 0
  stop = FALSE
  for (i in 1:NROW(F_M)){
    for (j in 1:NROW(Lc)){
      tmpPars<-MyPars
      tmpPars@FM<-F_M[i]
      tmpPars@SL50 <- Lc[j]
      tmpPars@SL95 <-Lc[j]+1
      tmpSim <- show_condition(LBSPRsim(tmpPars, verbose=FALSE))

      if(is.na(tmpSim)[1]) {
        stop = TRUE
        break
      }

      SPR_EU[i,j]=tmpSim$SPR
      YPR_EU[i,j]=tmpSim$YPR
      Yield_EU[i,j]=tmpSim$Yield

      counter<-counter+1

      if(!is.null(progressName)){
        shinyWidgets::updateProgressBar(
          session = sessionName,
          id = progressName,
          value = round(counter/steps*100,0)
        )
      }
    }
    if(stop) break
  }

  Yield_EU<- Yield_EU/max(Yield_EU, na.rm=TRUE)


  return(new("LBSPRarray",
             LifeHistory = LifeHistory,
             sim=list(Lc = Lc, F_M = F_M, SPR_EU = SPR_EU, YPR_EU = YPR_EU, Yield_EU=Yield_EU, stop = stop)
        )
  )
}
