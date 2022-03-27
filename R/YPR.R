

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
#'may want to re-calculate this quantity as relative YPR for presentation purposes.
#'Required parameters of the LifeHistory object are: Linf, L50, L95, MK, Steep, LW_A, LW_B
#' @param LifeHistoryObj  A life history object.
#' @param binWidth LBSPR length bin width, default value is 1.
#' @param binMin LBSPR binMin, default is 0
#' @param LcStep Length step size in cm for sequence of length at vulnerability. Approx. knife edge vul. SL50 = Lc, SL95 = Lc + 1
#' @param F_MStep F/M ratio step size for sequence of F_M
#' @param waitName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @param hostName When used within a shiny app, this function can update a host from the waiter package. See example.
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
#' #' sim<-lbsprSimWrapper(lh)
#'
#' #################################################
#' #Use of hostess loading bar from waiter package
#' ################################################
#' \dontrun{
#' library(shiny)
#' library(waiter)
#' library(fishSimGTG)
#'
#' ui <- fluidPage(
#'  useWaiter(),
#'   useHostess(), # include dependencies
#'   actionButton("btn", "render"),
#')
#'
#' server <- function(input, output){
#'
#'   host <- Hostess$new()
#'   w <- Waiter$new(
#'     html = host$get_loader(
#'       preset = "bubble",
#'       text_color = "black",
#'       center_page = TRUE,
#'       class = "",
#'       min = 0,
#'       max = 100,
#'       svg = NULL,
#'       progress_type = "fill",
#'       fill_direction = c("btt", "ttb", "ltr", "rtl"),
#'       stroke_direction = c("normal", "reverse"),
#'       fill_color = NULL,
#'       stroke_color = "pink"
#'     ),
#'     color = transparent(alpha = 0.2),
#'     fadeout = TRUE
#'   )
#'
#'  w$show()
#'   lbsprSimWrapper(LifeHistoryObj = LifeHistoryExample, waitName=w, hostName=host)
#'   w$hide()
#'
#' }
#'
#'
#' shinyApp(ui, server)}


lbsprSimWrapper<-function(LifeHistoryObj, binWidth=1, binMin=0, LcStep = 1, F_MStep = 0.2, waitName=NULL, hostName=NULL){

  if(!is.numeric(binWidth) ||
     !is.numeric(binMin) ||
     !is.numeric(LcStep) ||
     !is.numeric(F_MStep) ||
     binWidth < 0 ||
     binMin < 0 ||
     LcStep < 0 ||
     F_MStep < 0 ||
     length(LifeHistoryObj@Linf) == 0 ||
     length(LifeHistoryObj@L50) == 0 ||
     length(LifeHistoryObj@L95) == 0 ||
     length(LifeHistoryObj@MK) == 0 ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@L50 < 0 ||
     LifeHistoryObj@MK < 0 ||
     LifeHistoryObj@L50 >= LifeHistoryObj@Linf ||
     LifeHistoryObj@L50 >= LifeHistoryObj@L95
  ) {
    return(new("LBSPRarray",
               LifeHistory = LifeHistoryObj,
               sim=list(stop = TRUE))
    )
  } else {

    #-----------------------------
    #Create life history pars list
    #-----------------------------

    MyPars <- new("LB_pars")
    MyPars@Linf <- LifeHistoryObj@Linf
    MyPars@L50 <- LifeHistoryObj@L50
    MyPars@L95 <- LifeHistoryObj@L95
    MyPars@MK <- LifeHistoryObj@MK
    MyPars@BinWidth <- binWidth
    MyPars@BinMin <- binMin
    if(length(LifeHistoryObj@Steep) > 0) MyPars@Steepness<-LifeHistoryObj@Steep
    if(length(LifeHistoryObj@L_units) > 0) MyPars@L_units <- LifeHistoryObj@L_units
    if(length(LifeHistoryObj@LW_A) > 0) MyPars@Walpha <- LifeHistoryObj@LW_A
    if(length(LifeHistoryObj@Walpha_units) > 0) MyPars@Walpha_units <-LifeHistoryObj@Walpha_units
    if(length(LifeHistoryObj@LW_B) > 0)  MyPars@Wbeta <- MyPars@FecB <- LifeHistoryObj@LW_B

    #Setup place holder values for these parameters, we will change these later
    MyPars@SL50 <- LifeHistoryObj@L50
    MyPars@SL95 <- LifeHistoryObj@L95
    MyPars@FM<-1

    #------------------
    #Eumetric analysis
    #------------------

    Lmax<-(1 - 0.01^(1/MyPars@MK)) * MyPars@Linf
    Lc<-seq(floor(0.1*Lmax),  floor(Lmax), LcStep)
    F_M<-round(seq(0, 4, F_MStep), 3)
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
    if(!is.null(hostName) & !is.null(waitName)){
      waitName$show()
    }
    for (i in 1:NROW(F_M)){
      for (j in 1:NROW(Lc)){
        tmpPars<-MyPars
        tmpPars@FM<-F_M[i]
        tmpPars@SL50 <- Lc[j]
        tmpPars@SL95 <- Lc[j]+1
        tmpSim <- show_condition(LBSPRsim(tmpPars, verbose=FALSE))



        if(is.na(tmpSim)[1]) {
          stop = TRUE
          break
        }

        SPR_EU[i,j]=tmpSim$SPR
        YPR_EU[i,j]=tmpSim$YPR
        Yield_EU[i,j]=tmpSim$Yield

        counter<-counter+1

        if(!is.null(hostName) & !is.null(waitName)){
          hostName$set(counter/steps*100)
        }
      }
      if(stop) break
    }

    if(!is.null(hostName) & !is.null(waitName)){
      waitName$hide()
    }
    Yield_EU<- Yield_EU/max(Yield_EU, na.rm=TRUE)

    return(new("LBSPRarray",
               LifeHistory = LifeHistoryObj,
               sim=list(Lc = Lc, F_M = F_M, SPR_EU = SPR_EU, YPR_EU = YPR_EU, Yield_EU=Yield_EU, LcStep = LcStep, F_MStep = F_MStep, stop = stop))
    )
  }
}


#----------------------
#LBSPR sim wrapper
#----------------------

#Roxygen header
#'LBSPR using 'absel' methodeumetric simulation
#'
#'For diagnostic purposes only. See lbsprSimWrapper for details.
#'
#' @param LifeHistoryObj  A life history object.
#' @param binWidth LBSPR length bin width, default value is 1.
#' @param binMin LBSPR binMin, default is 0
#' @param LcStep Length step size in cm for sequence of length at vulnerability. Approx. knife edge vul. SL50 = Lc, SL95 = Lc + 1
#' @param F_MStep F/M ratio step size for sequence of F_M
#' @param waitName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @param hostName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @import LBSPR
#' @importFrom shinyWidgets updateProgressBar
#' @importFrom methods new
#' @export

lbsprSimWrapperAbsel<-function(LifeHistoryObj, binWidth=1, binMin=0, LcStep = 1, F_MStep = 0.2, waitName=NULL, hostName=NULL){

  if(!is.numeric(binWidth) ||
     !is.numeric(binMin) ||
     !is.numeric(LcStep) ||
     !is.numeric(F_MStep) ||
     binWidth < 0 ||
     binMin < 0 ||
     LcStep < 0 ||
     F_MStep < 0 ||
     length(LifeHistoryObj@Linf) == 0 ||
     length(LifeHistoryObj@L50) == 0 ||
     length(LifeHistoryObj@L95) == 0 ||
     length(LifeHistoryObj@MK) == 0 ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@L50 < 0 ||
     LifeHistoryObj@MK < 0 ||
     LifeHistoryObj@L50 >= LifeHistoryObj@Linf ||
     LifeHistoryObj@L50 >= LifeHistoryObj@L95
  ) {
    return(new("LBSPRarray",
               LifeHistory = LifeHistoryObj,
               sim=list(stop = TRUE))
    )
  } else {

    #-----------------------------
    #Create life history pars list
    #-----------------------------

    MyPars <- new("LB_pars")
    MyPars@Linf <- LifeHistoryObj@Linf
    MyPars@L50 <- LifeHistoryObj@L50
    MyPars@L95 <- LifeHistoryObj@L95
    MyPars@MK <- LifeHistoryObj@MK
    MyPars@BinWidth <- binWidth
    MyPars@BinMin <- binMin
    if(length(LifeHistoryObj@Steep) > 0) MyPars@Steepness<-LifeHistoryObj@Steep
    if(length(LifeHistoryObj@L_units) > 0) MyPars@L_units <- LifeHistoryObj@L_units
    if(length(LifeHistoryObj@LW_A) > 0) MyPars@Walpha <- LifeHistoryObj@LW_A
    if(length(LifeHistoryObj@Walpha_units) > 0) MyPars@Walpha_units <-LifeHistoryObj@Walpha_units
    if(length(LifeHistoryObj@LW_B) > 0)  MyPars@Wbeta <- MyPars@FecB <- LifeHistoryObj@LW_B

    #Setup place holder values for these parameters, we will change these later
    MyPars@SL50 <- LifeHistoryObj@L50
    MyPars@SL95 <- LifeHistoryObj@L95
    MyPars@FM<-1

    #------------------
    #Eumetric analysis
    #------------------

    Lmax<-(1 - 0.01^(1/MyPars@MK)) * MyPars@Linf
    Lc<-seq(floor(0.1*Lmax),  floor(Lmax), LcStep)
    F_M<-round(seq(0, 4, F_MStep), 3)
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
    if(!is.null(hostName) & !is.null(waitName)){
      waitName$show()
    }
    for (i in 1:NROW(F_M)){
      for (j in 1:NROW(Lc)){
        tmpPars<-MyPars
        tmpPars@FM<-F_M[i]
        tmpPars@SL50 <- Lc[j]
        tmpPars@SL95 <- Lc[j]+1
        tmpSim <- show_condition(LBSPRsim(tmpPars, verbose=FALSE, Control = list(modtype = "absel")))

        X<-LBSPRsim(tmpPars, verbose=FALSE, Control = list(modtype = "absel"))


        if(is.na(tmpSim)[1]) {
          stop = TRUE
          break
        }

        SPR_EU[i,j]=tmpSim$SPR
        YPR_EU[i,j]=tmpSim$YPR
        Yield_EU[i,j]=tmpSim$Yield

        counter<-counter+1

        if(!is.null(hostName) & !is.null(waitName)){
          hostName$set(counter/steps*100)
        }
      }
      if(stop) break
    }

    if(!is.null(hostName) & !is.null(waitName)){
      waitName$hide()
    }
    Yield_EU<- Yield_EU/max(Yield_EU, na.rm=TRUE)

    return(new("LBSPRarray",
               LifeHistory = LifeHistoryObj,
               sim=list(Lc = Lc, F_M = F_M, SPR_EU = SPR_EU, YPR_EU = YPR_EU, Yield_EU=Yield_EU, LcStep = LcStep, F_MStep = F_MStep, stop = stop))
    )
  }
}

#----------------------
#YPR fishSimGTG wrapper
#----------------------

#Roxygen header
#'YPR eumetric simulation using growth-type group model
#'
#'Produces YPR, Yield, and SPR arrays across combinations of F_M and Lc
#'
#'Utilizes the growth-type group model to calculate YPR, Yield, and SPR across factorial combination of F_M and Lc.
#'This produces an array of values for each output, which is useful for surface plot, finding MSY, etc.
#'Yield is presented in relative terms, with a maximum of 1. YPR is not presented in relative terms, but the user
#'may want to re-calculate this quantity as relative YPR for presentation purposes.
#'Required parameters of the LifeHistory object are: Linf, L50, L95, M, K, Steep, LW_A, LW_B
#' @param LifeHistoryObj  A life history object.
#' @param LcStep Length step size in cm for sequence of length at vulnerability. Approx. knife edge vul. SL50 = Lc, SL95 = Lc + 1
#' @param F_MStep F/M ratio step size for sequence of F_M
#' @param waitName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @param hostName When used within a shiny app, this function can update a host from the waiter package. See example.
#' @param gtg The number of growth-type groups. Default is 13.
#' @param stepsPerYear The number of steps per year. Default is 12.
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
#' #' sim<-gtgYPRWrapper(lh)
#'
#' #################################################
#' #Use of hostess loading bar from waiter package
#' ################################################
#' \dontrun{
#' library(shiny)
#' library(waiter)
#' library(fishSimGTG)
#'
#' ui <- fluidPage(
#'  useWaiter(),
#'   useHostess(), # include dependencies
#'   actionButton("btn", "render"),
#')
#'
#' server <- function(input, output){
#'
#'   host <- Hostess$new()
#'   w <- Waiter$new(
#'     html = host$get_loader(
#'       preset = "bubble",
#'       text_color = "black",
#'       center_page = TRUE,
#'       class = "",
#'       min = 0,
#'       max = 100,
#'       svg = NULL,
#'       progress_type = "fill",
#'       fill_direction = c("btt", "ttb", "ltr", "rtl"),
#'       stroke_direction = c("normal", "reverse"),
#'       fill_color = NULL,
#'       stroke_color = "pink"
#'     ),
#'     color = transparent(alpha = 0.2),
#'     fadeout = TRUE
#'   )
#'
#'  w$show()
#'   gtgYPRWrapper(LifeHistoryObj = LifeHistoryExample, waitName=w, hostName=host)
#'   w$hide()
#'
#' }
#'
#'
#' shinyApp(ui, server)}


gtgYPRWrapper<-function(LifeHistoryObj, LcStep = 1, F_MStep = 0.2, waitName=NULL, hostName=NULL, gtg=13, stepsPerYear=12){

  #-----------------------------
  #Initial check of conditions
  #-----------------------------
  if(!is.numeric(LcStep) ||
     !is.numeric(F_MStep) ||
     LcStep < 0 ||
     F_MStep < 0 ||
     class(LifeHistoryObj) != "LifeHistory"
  ) {
    return(new("YPRarray",
               lh = NULL,
               sim=list(stop = TRUE))
    )
  } else {

    #----------------------------------------------------------
    #LH assumptions and check of sufficient life history
    #If certain params not provided, then assumptions are made.
    #Same assumptions are made as in LBSPRsim for compatability
    #-----------------------------------------------------------
    if(length(LifeHistoryObj@Steep) == 0) LifeHistoryObj@Steep<-0.7
    if(LifeHistoryObj@Steep < 0.2) LifeHistoryObj@Steep<-0.7
    if(LifeHistoryObj@Steep > 1.0) LifeHistoryObj@Steep<-0.7
    if(length(LifeHistoryObj@LW_A) == 0) LifeHistoryObj@LW_A<-1e-04
    if(length(LifeHistoryObj@LW_B) == 0) LifeHistoryObj@LW_B<-3
    if(length(LifeHistoryObj@R0) == 0) LifeHistoryObj@R0<-10000

    FisheryObj<-new("Fishery")
    FisheryObj@vulType<-"logistic"
    FisheryObj@retType<-"full"
    FisheryObj@retMax<-1.0
    FisheryObj@Dmort<-0.0

    TimeAreaObj<-new("TimeArea")
    TimeAreaObj@gtg<-gtg

    #---------------------------------------------------
    #Check that life history requirements have been met
    #---------------------------------------------------
    lh<-LHwrapper(LifeHistoryObj, TimeAreaObj, stepsPerYear)
    if(is.null(lh)) {
      return(new("YPRarray",
                 lh = lh,
                 sim=list(stop = TRUE))
      )
    } else {

      #------------------
      #Eumetric analysis
      #------------------

      Lmax<-(1 - 0.01^(1/(LifeHistoryObj@M/LifeHistoryObj@K))) * LifeHistoryObj@Linf
      Lc<-seq(floor(0.1*Lmax),  floor(Lmax), LcStep)
      F_M<-round(seq(0, 4, F_MStep), 3)
      SPR_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))
      YPR_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))
      Yield_EU<-matrix(nrow=NROW(F_M), ncol=NROW(Lc))

      show_condition <- function(code) {
        tryCatch({
          code
        },
        error = function(c) NULL
        )
      }

      steps<-NROW(Lc)*NROW(F_M)
      counter <- 0
      stop = FALSE
      if(!is.null(hostName) & !is.null(waitName)){
        waitName$show()
      }
      for (j in 1:NROW(Lc)){

        FisheryObj@vulParams<-c(Lc[j], Lc[j]+1)
        sel<-selWrapper(lh, TimeAreaObj, FisheryObj, doPlot = FALSE)

        for (i in 1:NROW(F_M)){

          Feq<-F_M[i]*lh$LifeHistory@M
          tmpSim<-show_condition(solveD(lh = lh, sel=sel, doFit = FALSE, F_in = Feq))

          if(is.null(tmpSim)[1]) {
            stop = TRUE
            break
          }

          SPR_EU[i,j]=tmpSim$SPR
          YPR_EU[i,j]=tmpSim$YPR
          Yield_EU[i,j]=tmpSim$catchB

          counter<-counter+1

          if(!is.null(hostName) & !is.null(waitName)){
            hostName$set(counter/steps*100)
          }

        }
        if(stop) break
      }

      if(!is.null(hostName) & !is.null(waitName)){
        waitName$hide()
      }
      Yield_EU<- Yield_EU/max(Yield_EU, na.rm=TRUE)

      return(new("YPRarray",
                 lh = lh,
                 sim=list(Lc = Lc, F_M = F_M, SPR_EU = SPR_EU, YPR_EU = YPR_EU, Yield_EU=Yield_EU, LcStep = LcStep, F_MStep = F_MStep, stop = stop))
      )
    }
  }
}
