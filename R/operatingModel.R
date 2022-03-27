

#---------------------------------------
#Life history wrapper for sub-cohorts
#---------------------------------------

#Roxygen header
#'Life history wrapper for sub-cohorts
#'
#'Creates the necessary age-based vectors describing life history of sub-cohorts
#'
#' @param LifeHistoryObj  A life history object.
#' @param TimeAreaObj A time-area object
#' @importFrom methods slot slotNames
#' @importFrom stats dnorm plogis
#' @export
#' @examples
#' ta<-new("TimeArea")
#' ta@gtg<-13
#' ta@stepsPerYear<-1
#'LHwrapper(LifeHistoryObj = LifeHistoryExample, TimeAreaObj=ta)

LHwrapper<-function(LifeHistoryObj, TimeAreaObj, stepsPerYear = 1){

  if(class(LifeHistoryObj) != "LifeHistory" ||
     length(LifeHistoryObj@Linf) == 0 ||
     length(LifeHistoryObj@L50) == 0 ||
     length(LifeHistoryObj@L95) == 0 ||
     length(LifeHistoryObj@K) == 0 ||
     length(LifeHistoryObj@M) == 0 ||
     length(LifeHistoryObj@LW_A) == 0 ||
     length(LifeHistoryObj@LW_B) == 0 ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@L50 < 0 ||
     LifeHistoryObj@M < 0 ||
     LifeHistoryObj@K < 0 ||
     LifeHistoryObj@L50 >= LifeHistoryObj@Linf ||
     LifeHistoryObj@L50 >= LifeHistoryObj@L95 ||
     class(TimeAreaObj) != "TimeArea" ||
     length(TimeAreaObj@gtg) == 0 ||
     stepsPerYear < 1
  ) {
    NULL
  } else {

    #----------------
    #How many gtg?
    #----------------
    gtg<-ceiling(ifelse(TimeAreaObj@gtg < 3, 3, TimeAreaObj@gtg))
    gtg<-ifelse((gtg %% 2) == 0, gtg+1, gtg)
    CVLinf<-0.1
    maxsd<-2 #number of standard deviations from mean Linf
    SDLinf<-LifeHistoryObj@Linf*CVLinf
    gtg_Linf <- seq(from = LifeHistoryObj@Linf - maxsd * SDLinf, to = LifeHistoryObj@Linf + maxsd * SDLinf, length.out = gtg)

    #---------------
    #Rec proportions
    #---------------
    recProb<-dnorm(gtg_Linf, mean = LifeHistoryObj@Linf, sd = SDLinf)
    recProb<-recProb/sum(recProb)

    #---------------
    #Age classes
    #---------------
    stepsPerYear<-floor(stepsPerYear)
    if(length(LifeHistoryObj@Tmax) == 0 || LifeHistoryObj@Tmax < 2) {
      ages<-seq(0, ceiling(-log(0.01)/LifeHistoryObj@M),1/stepsPerYear)
    } else {
      ages<-seq(0, LifeHistoryObj@Tmax,1/stepsPerYear)
    }

    #---------------------
    #Life history vectors
    #---------------------
    t0<-ifelse(length(LifeHistoryObj@t0) == 0, 0, LifeHistoryObj@t0)
    L<-lapply(gtg_Linf, FUN=function(x) {
      tmp<-x*(1-exp(-LifeHistoryObj@K*(ages-t0)))
      tmp[tmp < 0] <- 0
      tmp
    })
    W<-lapply(1:gtg, FUN=function(x) LifeHistoryObj@LW_A*L[[x]]^LifeHistoryObj@LW_B)

    mat<-list()
    for (l in 1:gtg){
      if(length(LifeHistoryObj@isHermaph) != 0 && LifeHistoryObj@isHermaph){
        s<--(LifeHistoryObj@H95-LifeHistoryObj@H50)/log(1/0.95-1)
        probFemale<-(1-plogis(L[[l]], location=LifeHistoryObj@H50, scale=s))
      } else {
        probFemale<-0.5
      }
      s<--(LifeHistoryObj@L95-LifeHistoryObj@L50)/log(1/0.95-1)
      mat[[l]]<-plogis(L[[l]], location=LifeHistoryObj@L50, scale=s)*probFemale
    }

    return(list(
      LifeHistory = LifeHistoryObj,
      W=W,
      mat=mat,
      L=L,
      gtg = gtg,
      stepsPerYear = stepsPerYear,
      recProb = recProb,
      ageClasses = NROW(ages),
      maxAge=ages[NROW(ages)]))
  }
}


#---------------------------------------
#Selectivity wrapper for sub-cohorts
#---------------------------------------

#Roxygen header
#'Selectivity wrapper for sub-cohorts
#'
#'Creates the necessary age-based vectors describing selectivity of sub-cohorts
#'
#'Selectivity or vulnerability, retention and discards are described as follows.
#'Vulnerability is the probability of being selected or vulnerable to the fishing gear. Shape with respect to length, with a maximum value of 1 is enforced.
#'
#'Retention is the probability of being retained and not discarded. Describes both the shape with respect to length and maximum value which can be less than 1.
#'
#'Keep is resulting probability of being landed: Vul x Ret. Used in reporting the catch and catch-at-age
#'
#'Dead discards is: Vul x (1-Ret) x D, where D is the discard death rate. Used in reporting dead discards and dead discards-at-age
#'
#'Total removals is: vul x (Ret + (1 - Ret) x D)
#'
#'Selectivity types: "logistic" with params vector c(length at 50% sel, length at 95% sel)
#'
#'Retention types: "full" with no params, assumes Keep = Ret, no discards, no discard mortality; "logistic" with params vector c(length at 50% ret, length at 95% ret); "slotLimit" with params vector c(min length, max length) where catches occur betweem min and max
#'
#'Total dead is: Vul x (Ret + (1-Ret)D)
#' @param lh  An object produced by LHWrapper.
#' @param TimeAreaObj A time-area object
#' @param FisheryObj A stock object
#' @param doPlot Creates a basic plot to visualize outcomes. Useful for ensuring parameter selections are sensible.
#' @importFrom methods slot slotNames
#' @importFrom graphics par lines legend
#' @importFrom stats median
#' @export

selWrapper<-function(lh, TimeAreaObj, FisheryObj, doPlot = FALSE){

  #logistic
  logisticProb<-function(L, param, maxProb){
    if(
      length(param) != 2 ||
      param[1] < 0 ||
      param[1] >= param[2] ||
      length(maxProb) == 0 ||
      maxProb < 0 ||
      maxProb > 1
    ) {
     NULL
    } else {
      tryCatch({
        plogis(L, location=param[1], scale= -(param[2]-param[1])/log(1/0.95-1))*maxProb
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  #Full prob
  fullProb<-function(L, maxProb){
    rep(1.0, NROW(L))*maxProb
  }

  #Slot limit
  slotProb<-function(L, param, maxProb){
    if(
      length(param) != 2 ||
      param[1] < 0 ||
      param[1] >= param[2] ||
      length(maxProb) == 0 ||
      maxProb < 0 ||
      maxProb > 1
    ) {
      NULL
    } else {
      tryCatch({
        ifelse(L >= param[1] & L <= param[2], 1.0, 0)*maxProb
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  sel<-list()
  if(is.null(lh) ||
     class(FisheryObj) != "Fishery" ||
     !(FisheryObj@vulType %in%  "logistic") ||
     !(FisheryObj@retType %in%  c("full", "logistic", "slotLimit")) ||
     length(FisheryObj@retMax) == 0 ||
     FisheryObj@retMax < 0 ||
     FisheryObj@retMax > 1 ||
     length(FisheryObj@Dmort) == 0 ||
     FisheryObj@Dmort < 0 ||
     FisheryObj@Dmort > 1
  ) {
    sel<-NULL
  } else {
    #Vulnerability
    if(FisheryObj@vulType == "logistic") {
      sel$vul<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@vulParams, maxProb = 1.0))
    }

    #Retention
    if(FisheryObj@retType == "logistic") {
      historical$ret<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@retParams, maxProb = FisheryObj@retMax))
    }

    if(FisheryObj@retType == "full") {
      sel$ret<-lapply(1:lh$gtg, FUN=function(x) fullProb(L = lh$L[[x]], maxProb = FisheryObj@retMax))
    }

    if(FisheryObj@retType == "slotLimit") {
      sel$ret<-lapply(1:lh$gtg, FUN=function(x) slotProb(L = lh$L[[x]], param = FisheryObj@retParams, maxProb = FisheryObj@retMax))
    }

    #Do Keep, Dead discards, total removals
    if(
      list(NULL) %in% sel$vul ||
      list(NULL) %in% sel$ret
    ) {
      sel<-NULL
    } else {
      #Keep
      sel$keep<-lapply(1:lh$gtg, FUN=function(x) sel$vul[[x]] * sel$ret[[x]])

      #Dead discards
      sel$discard<-lapply(1:lh$gtg, FUN=function(x) sel$vul[[x]] * (1-sel$ret[[x]]) * FisheryObj@Dmort)

      #Removals (catch plus dead discards)
      sel$removal<-lapply(1:lh$gtg, FUN=function(x) sel$vul[[x]] * (sel$ret[[x]] + (1-sel$ret[[x]]) * FisheryObj@Dmort))
    }
  }

  #doPlot
  if(doPlot){
    if(!is.null(sel)){
      par(mfcol=c(1,1), las = 1)
      plot(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$vul)[order(unlist(lh$L))], type = "l", col = "purple", lwd =3, ylim = c(0,1), xlab = "Length", ylab = "Probability")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$ret)[order(unlist(lh$L))], lwd =3, col = "blue")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$keep)[order(unlist(lh$L))], lwd =3, col = "green", lty = 3)
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$discard)[order(unlist(lh$L))], lwd =3, col = "red")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$removal)[order(unlist(lh$L))], lwd =3, col = "orange", lty = 2)
      legend("topleft", legend = c("Vulnerability", "Retention", "Keep", "Dead discards", "Removals"), fill = c("purple", "blue", "green", "red", "orange"), border = "grey", bty = "n", inset=c(0, 0.1), cex = 0.8, x.intersp = 0.3)
    }
  }
  return(sel)
}


#-----------------------------------------
#Equilibrium calculations for sub-cohorts
#-----------------------------------------

#Roxygen header
#'Equilibrium conditions for sub-cohorts
#'
#'Creates the necessary age-based vectors equilibrium abundance, biomass and catch for sub-cohorts
#' @param lh  An object produced by LHWrapper.
#' @param sel An object produced by selWrapper
#' @param doFit Logical. When TRUE, estimates equilibrium fishing mortality based on input D_in. Ignores F_in. Default is FALSE
#' @param F_in Equilibrium fishing mortality rate. Used to calculate equilibrium conditions of the stock. Ignored when doFit = TRUE
#' @param D_type When doFit = TRUE, specifies type of equilibrium state metric that is specified in D_in (e.g., SSB depletion or SPR). Currently not in use and only SSB depletion is supported
#' @param D_in When doFit = TRUE, specifies value of equilibrium state. Currently this must be SSB depletion (value between 0 and 1)
#' @param doPlot Equilibrium length composition
#' @importFrom methods slot slotNames
#' @import ggplot2 gridExtra dplyr
#' @importFrom stats optimize
#' @export
#' @examples
#' sim<-gtgYPRWrapper(LifeHistoryObj = LifeHistoryExample, gtg=21, stepsPerYear = 12)

solveD<-function(lh, sel, doFit = FALSE, F_in = NULL, D_type = NULL, D_in = NULL, doPlot = FALSE){


  if(is.null(lh) ||
     is.null(sel) ||
     length(lh$LifeHistory@Steep) == 0 ||
     lh$LifeHistory@Steep < 0.21 ||
     lh$LifeHistory@Steep > 1 ||
     isTRUE(doFit & !(D_type %in%  c("relB", "SPR"))) ||
     isTRUE(doFit & is.null(D_in)) ||
     isTRUE(doFit & D_in < 0) ||
     isTRUE(doFit & D_in > 1)
  ) {
    return(NULL)
  } else {

    #---------------
    #Steps per year
    #---------------
    stepsPerYear <- lh$stepsPerYear
    totalSteps <- NROW(lh$L[[1]])

    #----------------------------------------
    #Fitting functions
    #----------------------------------------

    #FIX for loop so it runs fast
    #FIX fact that this function assumes age 1 rec, fix to reflect age 0 rec
    #FIX to include plus group

    min.Depletion<-function(logFmort){
      Fmort<-exp(logFmort)
      N<-lapply(1:lh$gtg, FUN=function(x) {
        tmp<-dplyr::lag(cumprod(exp(-lh$LifeHistory@M/stepsPerYear - Fmort/stepsPerYear*sel$removal[[x]])), n=1, default = 1)*lh$recProb[x]
        tmp[totalSteps]<- tmp[totalSteps]/(1-exp(-lh$LifeHistory@M/stepsPerYear - Fmort/stepsPerYear*sel$removal[[x]][totalSteps]))
        tmp
      })
      SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]]*lh$mat[[x]]*lh$W[[x]])[2:totalSteps])))
      SPR<- SB / Wbar
      D<-(4*lh$LifeHistory@Steep*SPR+lh$LifeHistory@Steep-1)/(5*lh$LifeHistory@Steep-1)
      if(D_type == "relB") return((D-D_in)^2)
      if(D_type == "SPR") return((SPR-D_in)^2)
    }

    #---------------
    #Wbar
    #---------------
    N<-lapply(1:lh$gtg, FUN=function(x) {
      tmp<-dplyr::lag(cumprod(rep(exp(-lh$LifeHistory@M/stepsPerYear), totalSteps)), n=1, default = 1)*lh$recProb[x]
      tmp[totalSteps]<- tmp[totalSteps]/(1-exp(-lh$LifeHistory@M/stepsPerYear))
      tmp
    })
    Wbar<-sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]]*lh$mat[[x]]*lh$W[[x]])[2:totalSteps])))

    #-------------
    #Get Feq
    #-------------
    if(doFit){
      Feq<-exp(optimize(min.Depletion, lower = -14, upper = 1.1,  maximum=FALSE, tol=0.00000001)$minimum)
    } else {
      Feq<-F_in
    }

    #-------------------------------------------------------------------------------------------
    #Get current SSB depletion
    #Re-cacl regardless of having D_in because in some cases fit cannot deplete stock to target
    #-------------------------------------------------------------------------------------------
    N<-lapply(1:lh$gtg, FUN=function(x) {
      tmp<-dplyr::lag(cumprod(exp(-lh$LifeHistory@M/stepsPerYear - Feq/stepsPerYear*sel$removal[[x]])), n=1, default = 1)*lh$recProb[x]
      tmp[totalSteps]<- tmp[totalSteps]/(1-exp(-lh$LifeHistory@M/stepsPerYear - Feq/stepsPerYear*sel$removal[[x]][totalSteps]))
      tmp
    })
    YPR<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq/stepsPerYear*sel$keep[[x]]/(Feq/stepsPerYear*sel$removal[[x]] + lh$LifeHistory@M/stepsPerYear)*(1-exp(-Feq/stepsPerYear*sel$removal[[x]]-lh$LifeHistory@M/stepsPerYear))*N[[x]])))
    SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]]*lh$mat[[x]]*lh$W[[x]])[2:totalSteps])))
    SPR<-SB / Wbar
    D<-max(0, (4*lh$LifeHistory@Steep*SPR+lh$LifeHistory@Steep-1)/(5*lh$LifeHistory@Steep-1))

    #------------------------------
    #Scale stock size and catches
    #------------------------------
    Req<-recruit(LifeHistoryObj = lh$LifeHistory, B0=Wbar, stock=Wbar*D, forceR=FALSE)
    N<-lapply(1:lh$gtg, FUN=function(x) {
      tmp<-dplyr::lag(cumprod(exp(-lh$LifeHistory@M/stepsPerYear - Feq/stepsPerYear*sel$removal[[x]])), n=1, default = 1)*lh$recProb[x]*Req
      tmp[totalSteps]<- tmp[totalSteps]/(1-exp(-lh$LifeHistory@M/stepsPerYear - Feq/stepsPerYear*sel$removal[[x]][totalSteps]))
      tmp
    })
    SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]]*lh$mat[[x]]*lh$W[[x]])[2:totalSteps])))
    VB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*sel$vul[[x]]*lh$W[[x]])))
    catchN<-sum(sapply(1:lh$gtg, FUN=function(x) sum(Feq*sel$keep[[x]]/(Feq*sel$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq/stepsPerYear*sel$removal[[x]]-lh$LifeHistory@M/stepsPerYear))*N[[x]])))
    catchB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq*sel$keep[[x]]/(Feq*sel$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq/stepsPerYear*sel$removal[[x]]-lh$LifeHistory@M/stepsPerYear))*N[[x]])))
    discN<-sum(sapply(1:lh$gtg, FUN=function(x) sum(Feq*sel$discard[[x]]/(Feq*sel$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq/stepsPerYear*sel$removal[[x]]-lh$LifeHistory@M/stepsPerYear))*N[[x]])))
    discB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq*sel$discard[[x]]/(Feq*sel$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq/stepsPerYear*sel$removal[[x]]-lh$LifeHistory@M/stepsPerYear))*N[[x]])))
    B0<-Wbar*lh$LifeHistory@R0

    if(doPlot) {
      tmp1<-data.frame(
        l = unlist(lapply(1:lh$gtg, FUN=function(x){
          rep(lh$L[[x]], N[[x]])
        })))
      tmp2<-data.frame(
        l = unlist(lapply(1:lh$gtg, FUN=function(x){
          rep(lh$L[[x]], N[[x]]*sel$vul[[x]])
        })))

      p1<-ggplot(tmp1, aes(x=l)) +
        geom_histogram(binwidth=1, position="identity", fill="#58C7B1", color="white") +
        labs(
          y= "Count",
          x = "Length",
          title = "Population length composition"
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(colour = "#808080"),
          legend.position = "right",
          strip.text.x = element_text(size=10, color = '#808080', face="bold"),
          strip.text.y = element_text(size=10, color = '#808080', face="bold"),
          axis.text.x = element_text( size = 10, color = '#808080', face="bold"),
          axis.text.y = element_text( size = 10, color = '#808080', face="bold"),
          axis.title.x = element_text(size=10, color = '#808080', face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size=10, color = '#808080', face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")),
          title = element_text(color = '#808080', face="bold")
        )

      p2<-ggplot(tmp2, aes(x=l)) +
        geom_histogram(binwidth=1, position="identity", fill="#58C7B1", color="white") +
        labs(
          y= "Count",
          x = "Length",
          title = "vulnerable length composition"
        ) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(colour = "#808080"),
          legend.position = "right",
          strip.text.x = element_text(size=10, color = '#808080', face="bold"),
          strip.text.y = element_text(size=10, color = '#808080', face="bold"),
          axis.text.x = element_text( size = 10, color = '#808080', face="bold"),
          axis.text.y = element_text( size = 10, color = '#808080', face="bold"),
          axis.title.x = element_text(size=10, color = '#808080', face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size=10, color = '#808080', face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")),
          title = element_text(color = '#808080', face="bold")
        )
      grid.arrange(p1, p2, nrow = 1)
    }

    #------------
    #Return list
    #------------
    return(list(Feq=Feq, D=D, SPR=SPR, Req=Req, B0=B0, SB=SB, VB=VB, catchN=catchN, catchB=catchB, discN=discN, discB=discB, N=N, YPR = YPR))
  }
}

#-----------------------------------------
#Stock - recruitment calculations
#-----------------------------------------

#Roxygen header
#'Stock-recuitment function
#'
#' @param LifeHistoryObj  A life history object.
#' @param B0 Unfished spawning stock biomass
#' @param stock the current spawning stock biomass
#' @param forceR A logical that overides the S-R function and simply returns a specified number of recruits. Default is FALSE
#' @param Rforced A fixed number of recruits to return. Only used when forceR is TRUE
#' @importFrom methods slot slotNames
#' @export

recruit<-function(LifeHistoryObj, B0=0, stock=0, forceR=FALSE, Rforced=0){
  R0<-LifeHistoryObj@R0
  h<-LifeHistoryObj@Steep
  if(forceR) Rnum<-Rforced
  if(!forceR) Rnum<-0.8*R0*h/(0.2*B0*(1-h)+(h-0.2)*stock)*stock
  return(Rnum)
}

#-----------------------------------------
#Rec deviations - recruitment calculations
#-----------------------------------------

#Roxygen header
#'Recuitment deviations
#'
#' @param LifeHistoryObj  A life history object.
#' @param TimeAreaObj A TimeArea object
#' @param StrategyObj A Strategy object
#' @importFrom methods slot slotNames
#' @export

recDev<-function(LifeHistoryObj, TimeAreaObj, StrategyObj = NULL){
  if(length(LifeHistoryObj@recSD) == 0 ||
     length(LifeHistoryObj@recRho) == 0 ||
     length(TimeAreaObj@historicalYears) == 0 ||
     length(TimeAreaObj@iterations) == 0 ||
     LifeHistoryObj@recSD < 0 ||
     LifeHistoryObj@recRho < 0 ||
     isTRUE(TimeAreaObj@historicalYears + ifelse(class(StrategyObj) == "Strategy" && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0) < 1) ||
     TimeAreaObj@iterations < 1
  ) {
    return(NULL)
  } else {
    years <- 1 + TimeAreaObj@historicalYears + ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
    iterations <- floor(TimeAreaObj@iterations)
    recSD <- LifeHistoryObj@recSD
    recRho <- LifeHistoryObj@recRho
    Rmult<-array(1:1, dim=c(years, iterations))
    for (k in 1:iterations){
      eps<-w<-rnorm(years,0,recSD)
      for (i in 2:NROW(eps)){
        eps[i]<-recRho*eps[i-1]+w[i]*sqrt(1-recRho*recRho)
      }
      Rmult[,k]<-exp(eps-recSD*recSD/2)
    }
    return(list(Rmult=Rmult))
  }
}


#-----------------------------------------
#Initial relative biomass deviations
#-----------------------------------------

#Roxygen header
#'Initial relative biomass deviations
#'
#' @param TimeAreaObj A TimeArea object
#' @param StochasticObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

bioDev<-function(TimeAreaObj, StochasticObj = NULL){
  if(length(TimeAreaObj@historicalBio) == 0 ||
     length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@historicalBio < 0 ||
     TimeAreaObj@historicalBio > 1 ||
     TimeAreaObj@iterations < 1
  ) {
    return(NULL)
  } else {

    if(class(StochasticObj) == "Stochastic" && length(StochasticObj@historicalBio) > 1 && StochasticObj@historicalBio[2] >= StochasticObj@historicalBio[1]) {
      Dtmp <- StochasticObj@historicalBio[1:2]
    } else {
      Dtmp <- c(TimeAreaObj@historicalBio, TimeAreaObj@historicalBio)
    }
    Dtmp<-ifelse(Dtmp < 0.01, 0.01, Dtmp) #Does your fishery even exist below 0.01?
    Dtmp<-ifelse(Dtmp > 0.95, 0.95, Dtmp) #Must start in fished state - Need initial F > 0
    iterations <- floor(TimeAreaObj@iterations)
    Ddev <- runif(n=iterations, min=Dtmp[1], max=Dtmp[2])
    return(list(Ddev=Ddev))
  }
}

#-----------------------------------------
#Initial cpue
#-----------------------------------------

#Roxygen header
#'Initial cpue deviations
#'
#' @param TimeAreaObj A TimeArea object
#' @param StochasticObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

cpueDev<-function(TimeAreaObj, StochasticObj = NULL){
  if(length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1
   ) {
    return(NULL)
  } else {
    if(class(StochasticObj) == "Stochastic" &&
       length(StochasticObj@historicalCPUE) > 1 &&
       StochasticObj@historicalCPUE[1] > 0 &&
       StochasticObj@historicalCPUE[2] > 0 &&
       StochasticObj@historicalCPUE[2] >= StochasticObj@historicalCPUE[1]
    ) {
      Ctmp <- StochasticObj@historicalCPUE[1:2]
    } else {
      Ctmp <- c(1, 1)
    }
    iterations <- floor(TimeAreaObj@iterations)
    Cdev <- runif(iterations, min = Ctmp[1], max = Ctmp[2])
    return(list(Cdev=Cdev))
  }
}


#-----------------------------------------
#Surival matrix calculations for cohort
#-----------------------------------------

#Roxygen header
#'Surival matrix calculations for cohort
#'
#' @param ageClasses Number of age classes. Typically max age + 1, since recruitment occurs at age 0
#' @param M_in Natural mortality rate
#' @param F_in Fishing mortality rate
#' @param S_in Selectivity, specifically the removal vector from selWrapper
#' @export

SurvMat<-function(ageClasses, M_in, F_in, S_in){
  Sout<-matrix(0:0, nrow=ageClasses, ncol=ageClasses)
  for(x in 1:ageClasses){
    Sout[x,x]<-exp(-M_in-S_in[x]*F_in)
  }

  Gout<-matrix(0:0, nrow=ageClasses, ncol=ageClasses)
  for(j in 1:(ageClasses-1)){
    Gout[(j+1),j]<-1
  }
  S<-Gout%*%Sout

  #Plus group
  S[ageClasses, ageClasses] <- exp(-M_in-S_in[ageClasses]*F_in)
  return(S)
}




#-----------------------------------------
#Move matrix calculations for cohort
#-----------------------------------------

#Roxygen header
#'Move matrix calculations for cohort
#'
#' @param Surv_in The output of the SurvMat function
#' @param Move_in A matrix of migration rates of dimensions areas x areas. Obtained from TimeArea object.
#' @param area_in The area under evaluation (e.g., 1)
#' @export

MoveMat<-function(Surv_in, Move_in, area_in){
  Mout<-Surv_in*Move_in[area_in,1]
  for(i in 2:NCOL(Move_in)){
    tmp<-Surv_in*Move_in[area_in,i]
    Mout<-cbind(Mout,tmp)
  }
  return(Mout)
}
