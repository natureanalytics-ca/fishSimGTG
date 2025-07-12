

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
#' @param stepsPerYear Creates multiple time steps per year allowing finer scale length composition.
#' @param doPlot Creates a basic plot to visualize outcomes. Useful for ensuring parameter selections are sensible.
#' @param wd A working directly where the output of runProjection is saved
#' @param imageName Character. A name for the resulting plot(s)
#' @param dpi Resolution in dots per inch of the resulting saved chart.
#' @importFrom methods slot slotNames
#' @importFrom stats dnorm plogis
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' ta<-new("TimeArea")
#' ta@gtg<-13
#'LHwrapper(LifeHistoryObj = LifeHistoryExample, TimeAreaObj=ta)

LHwrapper<-function(LifeHistoryObj, TimeAreaObj, stepsPerYear = 1, doPlot = FALSE, wd = NULL, imageName = NULL, dpi = 300){

  if(!is(LifeHistoryObj, "LifeHistory") ||
     length(LifeHistoryObj@Linf) == 0 ||
     length(LifeHistoryObj@L50) == 0 ||
     length(LifeHistoryObj@L95delta) == 0 ||
     length(LifeHistoryObj@K) == 0 ||
     length(LifeHistoryObj@M) == 0 ||
     length(LifeHistoryObj@LW_A) == 0 ||
     length(LifeHistoryObj@LW_B) == 0 ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@L50 < 0 ||
     LifeHistoryObj@M < 0 ||
     LifeHistoryObj@K < 0 ||
     LifeHistoryObj@L50 >= LifeHistoryObj@Linf ||
     isFALSE(LifeHistoryObj@L95delta > 0) ||
     !is(TimeAreaObj, "TimeArea") ||
     length(TimeAreaObj@gtg) == 0 ||
     stepsPerYear < 1
  ) {
    NULL
  } else {

    #----------------
    #How many gtg?
    #----------------
    gtg <- TimeAreaObj@gtg
    CVLinf<-0.1
    maxsd<-2 #number of standard deviations from mean Linf
    SDLinf<-LifeHistoryObj@Linf*CVLinf
    if(gtg == 1){
      gtg_Linf <- LifeHistoryObj@Linf
    } else {
      #gtg<-ceiling(ifelse(TimeAreaObj@gtg < 7, 7, TimeAreaObj@gtg))
      gtg<-ifelse((gtg %% 2) == 0, gtg+1, gtg)
      gtg_Linf <- seq(from = LifeHistoryObj@Linf - maxsd * SDLinf, to = LifeHistoryObj@Linf + maxsd * SDLinf, length.out = gtg)
    }

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
        s<--(LifeHistoryObj@H95delta)/log(1/0.95-1)
        probFemale<-(1-plogis(L[[l]], location=LifeHistoryObj@H50, scale=s))
      } else {
        probFemale<-0.5
      }
      s<--(LifeHistoryObj@L95delta)/log(1/0.95-1)
      mat[[l]]<-plogis(L[[l]], location=LifeHistoryObj@L50, scale=s)*probFemale
    }

    #----
    #Plot
    #----
    if(doPlot){
      md<-ceiling(gtg/2)
      dt<-data.frame(
        ages = ages,
        mat = mat[[md]],
        L = L[[md]],
        W = W[[md]]
      )
      p1<- ggplot(dt, aes(x = ages, y = L)) +
        geom_line(color = "cornflowerblue", size = 1.5) +
        ylab("Length") +
        xlab("Age") +
        theme_classic() +
        theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
              strip.text.y=element_text(colour = "black", size=8, face="bold"),
              strip.background = element_rect(fill ="lightgrey"),
              axis.text=element_text(size=8),
              panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))
      p2<- ggplot(dt, aes(x = ages, y = mat)) +
        geom_line(color = "cornflowerblue", size = 1.5) +
        ylab("Proportion mature females") +
        xlab("Age") +
        theme_classic() +
        theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
              strip.text.y=element_text(colour = "black", size=8, face="bold"),
              strip.background = element_rect(fill ="lightgrey"),
              axis.text=element_text(size=8),
              panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))
      p3<- ggplot(dt, aes(x = ages, y = W)) +
        geom_line(color = "cornflowerblue", size = 1.5) +
        ylab("Weight") +
        xlab("Age") +
        theme_classic() +
        theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
              strip.text.y=element_text(colour = "black", size=8, face="bold"),
              strip.background = element_rect(fill ="lightgrey"),
              axis.text=element_text(size=8),
              panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))
      p4<- ggplot(dt, aes(x = L, y = W)) +
        geom_line(color = "cornflowerblue", size = 1.5) +
        ylab("Weight") +
        xlab("Length") +
        theme_classic() +
        theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
              strip.text.y=element_text(colour = "black", size=8, face="bold"),
              strip.background = element_rect(fill ="lightgrey"),
              axis.text=element_text(size=8),
              panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))
      if(is.null(wd) | is.null(imageName)) grid.arrange(p1, p2, p3, p4, ncol = 2)
      if(!is.null(wd) & !is.null(imageName)) {
        p5<-grid.arrange(p1, p2, p3, p4, ncol = 2)
        ggsave(filename = paste0(wd, "/", imageName, "_LH.png"), plot = p5, device = "png", dpi = dpi, width = 6, height = 6, units = "in")
      }
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
#'Selectivity types: "logistic" with params vector c(length at 50% sel, length increment to 95% sel); "explog" exponential logistic (domed) with vector c(p1, peak, p3); "gillnetMasterNormal" master curve for mesh sel with vector(R0 relative peak L/m, var, mesh (in cm)); "gillnetMasterLognormal" master curve for mesh sel with vector (R0 relative peak L/m, var, mesh (in cm))
#'
#'Retention types: "full" with no params, assumes Keep = Ret; "logistic" with params vector c(length at 50% ret, length increment to 95% ret); "slotLimit" with params vector c(min length, max length) where catches occur betweem min and max
#'
#'Total dead is: Vul x (Ret + (1-Ret)D)
#' @param lh  An object produced by LHWrapper.
#' @param TimeAreaObj A time-area object
#' @param FisheryObj A fishery object
#' @param doPlot Creates a basic plot to visualize outcomes. Useful for ensuring parameter selections are sensible.
#' @param wd A working directly where the output of runProjection is saved
#' @param imageName Character. A name for the resulting plot(s)
#' @param dpi Resolution in dots per inch of the resulting saved chart.
#' @importFrom methods slot slotNames
#' @importFrom graphics par lines legend
#' @importFrom stats median
#' @export

selWrapper<-function(lh, TimeAreaObj, FisheryObj, doPlot = FALSE,  wd = NULL, imageName = NULL, dpi = 300){

  #logistic
  logisticProb<-function(L, param, maxProb){
    if(
      length(param) != 2 ||
      param[1] < 0 ||
      param[2] < 0 ||
      length(maxProb) == 0 ||
      maxProb < 0 ||
      maxProb > 1
    ) {
     NULL
    } else {
      tryCatch({
        plogis(L, location=param[1], scale= -(param[2])/log(1/0.95-1))*maxProb
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  #Exponential logistic
  explogProb<-function(L, param, maxProb){
    if(
      length(param) != 3 ||
      param[1] < 0.02 ||
      param[1] > 1 ||
      param[3] < 0.001 ||
      param[3] > 0.5 ||
      length(maxProb) == 0 ||
      maxProb < 0 ||
      maxProb > 1
    ) {
      NULL
    } else {
      tryCatch({
        exp(param[3]*param[1]*(param[2]-L))/(1-param[3]*(1-exp(param[1]*(param[2]-L))))
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  #Gillnet master normal
  gillnetNormalProb<-function(L, param){
    if(
      length(param) != 3 ||
      param[1] < 0 ||
      param[2] < 0 ||
      param[3] < 0
    ) {
      NULL
    } else {
      tryCatch({
        R<-L/param[3]
        exp(-(R - param[1])^2/(2*param[2]))
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  #Gillnet master lognormal
  gillnetLogProb<-function(L, param){
    if(
      length(param) != 3 ||
      param[1] < 0 ||
      param[2] < 0 ||
      param[3] < 0
    ) {
      NULL
    } else {
      tryCatch({
        R<-L/param[3]
        exp(-(log(R) - log(param[1]))^2/(2*param[2]))
      },
      error = function(c) NULL,
      warning = function(c) NULL
      )
    }
  }

  #Gillnet master skewed normal
  gillnetLogProb<-function(L, param){
    if(
      length(param) != 3 ||
      param[1] < 0 ||
      param[2] < 0 ||
      param[3] < 0
    ) {
      NULL
    } else {
      tryCatch({
        R<-L/param[3]
        exp(-(log(R) - log(param[1]))^2/(2*param[2]))
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
     !is(FisheryObj, "Fishery") ||
     !(FisheryObj@vulType %in%  c("logistic", "explog", "gillnetMasterNormal", "gillnetMasterLognormal")) ||
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
    if(FisheryObj@vulType == "explog") {
      sel$vul<-lapply(1:lh$gtg, FUN=function(x) explogProb(L = lh$L[[x]], param = FisheryObj@vulParams, maxProb = 1.0))
    }
    if(FisheryObj@vulType == "gillnetMasterNormal") {
      sel$vul<-lapply(1:lh$gtg, FUN=function(x) gillnetNormalProb(L = lh$L[[x]], param = FisheryObj@vulParams))
    }
    if(FisheryObj@vulType == "gillnetMasterLognormal") {
      sel$vul<-lapply(1:lh$gtg, FUN=function(x) gillnetLogProb(L = lh$L[[x]], param = FisheryObj@vulParams))
    }


    #Retention
    if(FisheryObj@retType == "logistic") {
      sel$ret<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@retParams, maxProb = FisheryObj@retMax))
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

    if(is.null(wd) | is.null(imageName)){
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
    if(!is.null(wd) & !is.null(imageName)) {
      if(!is.null(sel)){
        png(filename=paste0(wd, "/", imageName, ".png"), width=6, height=6, units="in", res=dpi, bg="white", pointsize=12)
        par(mfcol=c(1,1), las = 1)
        plot(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$vul)[order(unlist(lh$L))], type = "l", col = "purple", lwd =3, ylim = c(0,1), xlab = "Length", ylab = "Probability")
        lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$ret)[order(unlist(lh$L))], lwd =3, col = "blue")
        lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$keep)[order(unlist(lh$L))], lwd =3, col = "green", lty = 3)
        lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$discard)[order(unlist(lh$L))], lwd =3, col = "red")
        lines(unlist(lh$L)[order(unlist(lh$L))], unlist(sel$removal)[order(unlist(lh$L))], lwd =3, col = "orange", lty = 2)
        legend("topleft", legend = c("Vulnerability", "Retention", "Keep", "Dead discards", "Removals"), fill = c("purple", "blue", "green", "red", "orange"), border = "grey", bty = "n", inset=c(0, 0.1), cex = 0.8, x.intersp = 0.3)
        dev.off()
      }
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
#' @param D_type When doFit = TRUE, specifies type of equilibrium state metric that is specified in D_in (e.g., SSB depletion or SPR).
#' @param D_in When doFit = TRUE, specifies value of equilibrium state. Must be SSB depletion or SPR both with value between 0 and 1
#' @param doPlot Equilibrium length composition
#' @importFrom methods slot slotNames
#' @import ggplot2  dplyr
#' @importFrom stats optimize
#' @importFrom gridExtra grid.arrange
#' @export

solveD<-function(lh, sel, doFit = FALSE, F_in = NULL, D_type = NULL, D_in = NULL, doPlot = FALSE){

  l <- NULL

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
#' @importFrom stats rnorm
#' @export

recDev<-function(LifeHistoryObj, TimeAreaObj, StochasticObj, StrategyObj = NULL){
  if(length(LifeHistoryObj@recSD) == 0 ||
     length(LifeHistoryObj@recRho) == 0 ||
     length(TimeAreaObj@historicalYears) == 0 ||
     length(TimeAreaObj@iterations) == 0 ||
     LifeHistoryObj@recSD < 0 ||
     LifeHistoryObj@recRho < 0 ||
     LifeHistoryObj@recRho > 1 ||
     isTRUE(TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy") && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0) < 1) ||
     TimeAreaObj@iterations < 1
  ) {
    return(NULL)
  } else {


    #--------
    #recSD
    #--------
    iterations <- floor(TimeAreaObj@iterations)
    recSD <- rep(LifeHistoryObj@recSD, iterations)
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@recSD) > 1 &&
       StochasticObj@recSD[1] > 0 &&
       StochasticObj@recSD[2] > 0 &&
       StochasticObj@recSD[2] >= StochasticObj@recSD[1]
    ) {
      recSD<-runif(iterations, min = StochasticObj@recSD[1], max = StochasticObj@recSD[2])
    }


    #--------
    #recRho
    #--------
    recRho <- rep(LifeHistoryObj@recRho, iterations)
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@recRho) > 1 &&
       StochasticObj@recRho[1] >= 0 &&
       StochasticObj@recRho[1] <= 1 &&
       StochasticObj@recRho[2] >= 0 &&
       StochasticObj@recRho[2] <= 1 &&
       StochasticObj@recRho[2] >= StochasticObj@recRho[1]
    ) {
      recRho<-runif(iterations, min = StochasticObj@recRho[1], max = StochasticObj@recRho[2])
    }


    years <- 1 + TimeAreaObj@historicalYears + ifelse(is(StrategyObj, "Strategy")  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
    Rmult<-array(1:1, dim=c(years, iterations))
    for (k in 1:iterations){
      eps<-w<-rnorm(years,0,recSD[k])
      for (i in 2:NROW(eps)){
        eps[i]<-recRho[k]*eps[i-1]+w[i]*sqrt(1-recRho[k]*recRho[k])
      }
      Rmult[,k]<-exp(eps-recSD[k]*recSD[k]/2)
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
#' @importFrom stats runif
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

    if(is(StochasticObj, "Stochastic") && length(StochasticObj@historicalBio) > 1 && StochasticObj@historicalBio[2] >= StochasticObj@historicalBio[1]) {
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
#Initial cpue - special function for projectionStrategy
#-----------------------------------------

#Roxygen header
#'Initial cpue deviations - special function for projectionStrategy
#'
#' @param TimeAreaObj A TimeArea object
#' @param StrategyObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

cpueDev<-function(TimeAreaObj, StrategyObj){
  if(length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1
   ) {
    return(NULL)
  } else {
    if(is(StrategyObj, "Strategy") &&
       length(StrategyObj@projectionParams[['CPUE']]) > 1 &&
       StrategyObj@projectionParams[['CPUE']][1] > 0 &&
       StrategyObj@projectionParams[['CPUE']][2] > 0 &&
       StrategyObj@projectionParams[['CPUE']][2] >= StrategyObj@projectionParams[['CPUE']][1]
    ) {
      Ctmp <- StrategyObj@projectionParams[['CPUE']][1:2]
      iterations <- floor(TimeAreaObj@iterations)
      Cdev <- runif(iterations, min = Ctmp[1], max = Ctmp[2])
      return(list(Cdev=Cdev))
    } else {
      return(NULL)
    }
  }
}


#--------------------------------------------------------------------------------------
#Effort implementation error in projection - special function for projectionStrategy
#--------------------------------------------------------------------------------------

#Roxygen header
#'Effort implementation error in projection - special function for projectionStrategy
#'
#' @param TimeAreaObj A TimeArea object
#' @param StrategyObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

effortImpErrorDev<-function(TimeAreaObj, StrategyObj){
  if(length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1
  ) {
    return(NULL)
  } else {
    iterations <- floor(TimeAreaObj@iterations)
    if(is(StrategyObj, "Strategy") &&
       length(StrategyObj@projectionParams[['effortImpError']]) == 2 &&
       StrategyObj@projectionParams[['effortImpError']][1] >= 0 &&
       StrategyObj@projectionParams[['effortImpError']][2] >= 0 &&
       StrategyObj@projectionParams[['effortImpError']][2] >= StrategyObj@projectionParams[['effortImpError']][1]
    ) {
      Edev <- runif(iterations, min = StrategyObj@projectionParams[['effortImpError']][1], max = StrategyObj@projectionParams[['effortImpError']][2])
      return(list(Edev=Edev))
    } else {
      return(list(Edev=runif(iterations, min = 1, max =1)))
    }
  }
}


#-----------------------------------------
#Set of life history params
#-----------------------------------------

#Roxygen header
#'Set of life history params
#'
#' @param TimeAreaObj A TimeArea object
#' @param StochasticObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

lifehistoryDev<-function(TimeAreaObj, StochasticObj){
  if(length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1
  ) {
    return(NULL)
  } else {

    #Any life history parameter can have uncertainty.
    #For those parameters specified in Stochastic object, uniform sampling occurs (1 draw per iteration)
    #Otherwise, constants are obtained from LifeHistory object.

    iterations <- floor(TimeAreaObj@iterations)

    #--------
    #Linf
    #--------
    Linf<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@Linf) > 1 &&
       StochasticObj@Linf[1] > 0 &&
       StochasticObj@Linf[2] > 0 &&
       StochasticObj@Linf[2] >= StochasticObj@Linf[1]
    ) {
      Linf<-runif(iterations, min = StochasticObj@Linf[1], max = StochasticObj@Linf[2])
    }

    #--------
    #K
    #--------
    K<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@K) > 1 &&
       StochasticObj@K[1] > 0 &&
       StochasticObj@K[2] > 0 &&
       StochasticObj@K[2] >= StochasticObj@K[1]
    ) {
      K<-runif(iterations, min = StochasticObj@K[1], max = StochasticObj@K[2])
    }

    #--------
    #L50
    #--------
    L50<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@L50) > 1 &&
       StochasticObj@L50[1] > 0 &&
       StochasticObj@L50[2] > 0 &&
       StochasticObj@L50[2] >= StochasticObj@L50[1]
    ) {
      L50<-runif(iterations, min = StochasticObj@L50[1], max = StochasticObj@L50[2])
    }

    #--------
    #L95delta
    #--------
    L95delta<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@L95delta) > 1 &&
       StochasticObj@L95delta[1] > 0 &&
       StochasticObj@L95delta[2] > 0 &&
       StochasticObj@L95delta[2] >= StochasticObj@L95delta[1]
    ) {
      L95delta<-runif(iterations, min = StochasticObj@L95delta[1], max = StochasticObj@L95delta[2])
    }

    #--------
    #M
    #--------
    M<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@M) > 1 &&
       StochasticObj@M[1] > 0 &&
       StochasticObj@M[2] > 0 &&
       StochasticObj@M[2] >= StochasticObj@M[1]
    ) {
      M<-runif(iterations, min = StochasticObj@M[1], max = StochasticObj@M[2])
    }

    #--------
    #Steep
    #--------
    Steep<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@Steep) > 1 &&
       StochasticObj@Steep[1] >= 0.21 &&
       StochasticObj@Steep[1] <= 1 &&
       StochasticObj@Steep[2] >= 0.21 &&
       StochasticObj@Steep[2] <= 1 &&
       StochasticObj@Steep[2] >= StochasticObj@Steep[1]
    ) {
      Steep<-runif(iterations, min = StochasticObj@Steep[1], max = StochasticObj@Steep[2])
    }

    #--------
    #H50
    #--------
    H50<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@H50) > 1 &&
       StochasticObj@H50[1] > 0 &&
       StochasticObj@H50[2] > 0 &&
       StochasticObj@H50[2] >= StochasticObj@H50[1]
    ) {
      H50<-runif(iterations, min = StochasticObj@H50[1], max = StochasticObj@H50[2])
    }

    #--------
    #H95delta
    #--------
    H95delta<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@H95delta) > 1 &&
       StochasticObj@H95delta[1] > 0 &&
       StochasticObj@H95delta[2] > 0 &&
       StochasticObj@H95delta[2] >= StochasticObj@H95delta[1]
    ) {
      H95delta<-runif(iterations, min = StochasticObj@H95delta[1], max = StochasticObj@H95delta[2])
    }

    return(list(Linf=Linf, K=K, L50=L50, L95delta=L95delta, M=M, Steep=Steep, H50=H50, H95delta=H95delta))
  }
}

#-----------------------------------------
#Set of selectivity params
#-----------------------------------------

#Roxygen header
#'Set of selectivity params
#'
#' @param TimeAreaObj A TimeArea object
#' @param HistFisheryObj A fishery object for historical period
#' @param ProFisheryObj_list A fishery object list for projection period
#' @param StochasticObj A Stochastic object
#' @importFrom methods slot slotNames
#' @export

selDev<-function(TimeAreaObj, HistFisheryObj, ProFisheryObj_list=NULL, StochasticObj){
  if(length(TimeAreaObj@iterations) == 0 ||
     TimeAreaObj@iterations < 1 ||
     length(TimeAreaObj@areas) == 0 ||
     TimeAreaObj@areas < 2
  ) {
    return(NULL)
  } else {

    #Any life history parameter can have uncertainty.
    #For those parameters specified in Stochastic object, uniform sampling occurs (1 draw per iteration)
    #Otherwise, constants are obtained from LifeHistory object.

    iterations <- floor(TimeAreaObj@iterations)

    #------------------------------
    #Historical period vulnerability
    #------------------------------
    historical_vul<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@histFisheryVul) > 0 &&
       dim(StochasticObj@histFisheryVul)[1] == 2 &&
       sum(sapply(1:dim(StochasticObj@histFisheryVul)[2], function(x){StochasticObj@histFisheryVul[2,x] >= StochasticObj@histFisheryVul[1,x]})) == dim(StochasticObj@histFisheryVul)[2]
    ) {
      historical_vul<-sapply(1:dim(StochasticObj@histFisheryVul)[2], function(x){
        runif(iterations, min = StochasticObj@histFisheryVul[1,x], max = StochasticObj@histFisheryVul[2,x])
      })
    }

    #------------------------------
    #Projection period vulnerability
    #------------------------------
    projection_vul<-lapply(1:TimeAreaObj@areas, function(x){
      NULL
    })
    #Check to see if same selectivity elements should be applied to projection period
    if(
      is(StochasticObj, "Stochastic") &&
      length(StochasticObj@sameFisheryVul) > 0  &&
      StochasticObj@sameFisheryVul
    ) {
      projection_vul<-lapply(1:TimeAreaObj@areas, function(x){
        historical_vul
      })
    #Otherwise calculate unique projection selectivity elements
    } else {
      if(is(StochasticObj, "Stochastic") &&
         length(StochasticObj@proFisheryVul_list) > 0
      ) {
        for(i in 1:TimeAreaObj@areas){
          if(
            dim(StochasticObj@proFisheryVul_list[[i]])[1] == 2 &&
            sum(sapply(1:dim(StochasticObj@proFisheryVul_list[[i]])[2], function(x){StochasticObj@proFisheryVul_list[[i]][2,x] >= StochasticObj@proFisheryVul_list[[i]][1,x]})) == dim(StochasticObj@proFisheryVul_list[[i]])[2]
          ) {
            projection_vul[[i]]<-sapply(1:dim(StochasticObj@proFisheryVul_list[[i]])[2], function(x){
              runif(iterations, min = StochasticObj@proFisheryVul_list[[i]][1,x], max = StochasticObj@proFisheryVul_list[[i]][2,x])
            })
          }
        }
      }
    }

    #----------------------
    #Historical retention
    #----------------------
    histical_retention<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@histFisheryRet) > 0 &&
       dim(StochasticObj@histFisheryRet)[1] == 2 &&
       sum(sapply(1:dim(StochasticObj@histFisheryRet)[2], function(x){StochasticObj@histFisheryRet[2,x] >= StochasticObj@histFisheryRet[1,x]})) == dim(StochasticObj@histFisheryRet)[2]
    ) {
      histical_retention<-sapply(1:dim(StochasticObj@histFisheryRet)[2], function(x){
        runif(iterations, min = StochasticObj@histFisheryRet[1,x], max = StochasticObj@histFisheryRet[2,x])
      })
    }

    #------------------------------
    #Projection period retention
    #------------------------------
    projection_retention<-lapply(1:TimeAreaObj@areas, function(x){
      NULL
    })
    #Check to see if same retention elements should be applied to projection period
    if(
      is(StochasticObj, "Stochastic") &&
      length(StochasticObj@sameFisheryRet) > 0 &&
      StochasticObj@sameFisheryRet
    ) {
      projection_retention<-lapply(1:TimeAreaObj@areas, function(x){
        histical_retention
      })
      #Otherwise calculate unique projection retention elements
    } else {
      if(is(StochasticObj, "Stochastic") &&
         length(StochasticObj@proFisheryRet_list) > 0
      ) {
        for(i in 1:TimeAreaObj@areas){
          if(
            dim(StochasticObj@proFisheryRet_list[[i]])[1] == 2 &&
            sum(sapply(1:dim(StochasticObj@proFisheryRet_list[[i]])[2], function(x){StochasticObj@proFisheryRet_list[[i]][2,x] >= StochasticObj@proFisheryRet_list[[i]][1,x]})) == dim(StochasticObj@proFisheryRet_list[[i]])[2]
          ) {
            projection_retention[[i]]<-sapply(1:dim(StochasticObj@proFisheryRet_list[[i]])[2], function(x){
              runif(iterations, min = StochasticObj@proFisheryRet_list[[i]][1,x], max = StochasticObj@proFisheryRet_list[[i]][2,x])
            })
          }
        }
      }
    }

    #-----------------------
    #Historical Dmort
    #-----------------------
    historical_Dmort<-NULL
    if(is(StochasticObj, "Stochastic") &&
       length(StochasticObj@histFisheryDmort) > 0 &&
       dim(StochasticObj@histFisheryDmort)[1] == 2 &&
       sum(StochasticObj@histFisheryDmort <= 1) == (dim(StochasticObj@histFisheryDmort)[1]*dim(StochasticObj@histFisheryDmort)[2]) &&
       sum(StochasticObj@histFisheryDmort >= 0) == (dim(StochasticObj@histFisheryDmort)[1]*dim(StochasticObj@histFisheryDmort)[2]) &&
       sum(sapply(1:dim(StochasticObj@histFisheryDmort)[2], function(x){StochasticObj@histFisheryDmort[2,x] >= StochasticObj@histFisheryDmort[1,x]})) == dim(StochasticObj@histFisheryDmort)[2]
    ) {
      historical_Dmort<-sapply(1:dim(StochasticObj@histFisheryDmort)[2], function(x){
        runif(iterations, min = StochasticObj@histFisheryDmort[1,x], max = StochasticObj@histFisheryDmort[2,x])
      })
    }

    #-----------------------
    #Projection Dmort
    #-----------------------
    projection_Dmort<-lapply(1:TimeAreaObj@areas, function(x){
      NULL
    })
    if(
      is(StochasticObj, "Stochastic") &&
      length(StochasticObj@sameFisheryDmort) > 0 &&
      StochasticObj@sameFisheryDmort
    ) {
      projection_Dmort<-lapply(1:TimeAreaObj@areas, function(x){
        historical_Dmort
      })
      #Otherwise calculate unique projection retention elements
    } else {
      if(is(StochasticObj, "Stochastic") &&
         length(StochasticObj@proFisheryDmort_list) > 0
      ) {
        for(i in 1:TimeAreaObj@areas){
          if(
            dim(StochasticObj@proFisheryDmort_list[[i]])[1] == 2 &&
            sum(StochasticObj@proFisheryDmort_list[[i]] <= 1) == (dim(StochasticObj@proFisheryDmort_list[[i]])[1]*dim(StochasticObj@proFisheryDmort_list[[i]])[2]) &&
            sum(StochasticObj@proFisheryDmort_list[[i]] >= 0) == (dim(StochasticObj@proFisheryDmort_list[[i]])[1]*dim(StochasticObj@proFisheryDmort_list[[i]])[2]) &&
            sum(sapply(1:dim(StochasticObj@proFisheryDmort_list[[i]])[2], function(x){StochasticObj@proFisheryDmort_list[[i]][2,x] >= StochasticObj@proFisheryDmort_list[[i]][1,x]})) == dim(StochasticObj@proFisheryDmort_list[[i]])[2]
          ) {
            projection_Dmort[[i]]<-sapply(1:dim(StochasticObj@proFisheryDmort_list[[i]])[2], function(x){
              runif(iterations, min = StochasticObj@proFisheryDmort_list[[i]][1,x], max = StochasticObj@proFisheryDmort_list[[i]][2,x])
            })
          }
        }
      }
    }

    return(list(
      hist = list(vulParams = historical_vul, retParams = histical_retention, Dmort = historical_Dmort),
      pro = lapply(1:TimeAreaObj@areas, function(x){
        list(vulParams = projection_vul[[x]], retParams = projection_retention[[x]], Dmort = projection_Dmort[[x]])
      })
    ))
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



#-----------------------------------------
# Obs model indices
#-----------------------------------------

calculateIndex  <- function(simulation_result, IndexObj){

  # define the dimensions (to get the structure of the simulation)
  years <- dim(simulation_result$dynamics$VB)[1]       #total years hist+future
  iterations <- dim(simulation_result$dynamics$VB)[2]  #total iterations
  total_areas <- dim(simulation_result$dynamics$VB)[3] #total areas
  historicalYears <- simulation_result$TimeAreaObj@historicalYears #toal historical years- before management
  historical_end <- 1 + historicalYears

  # counts number of individual survey/CPUE programs defined in the design list.
  n_indices <- length(IndexObj@survey_design)

  if(n_indices == 0) {
    stop("survey_design list must contain at least one survey or CPUE design")
  }

  #  add validation for each survey/CPUE design
  #  lopp through each individual survey/CPUE design for validation.
  for(index_idx in 1:n_indices) {
    design <- IndexObj@survey_design[[index_idx]] #extracts current survey design from the list (for each index)

    # common required elements for surveys/CPUEs
    required_elements <- c("indextype", "areas", "indexYears",
                           "q_hist_bounds", "q_proj_bounds",
                           "hyperstability_hist_bounds", "hyperstability_proj_bounds",
                           "obsError_CV_hist_bounds", "obsError_CV_proj_bounds")

    # for FI surveys only, selectivity indices and survey_timing are needed
    if("indextype" %in% names(design) && design$indextype == "FI") {
      required_elements <- c(required_elements, "selectivity_hist_idx",
                             "selectivity_proj_idx", "survey_timing")
    }

    # check the required elements are available
    if(!all(required_elements %in% names(design))) { # check if required element exist
      missing <- required_elements[!required_elements %in% names(design)]   # identify what element is misssing
      stop(paste("Survey design", index_idx, "missing elements:", paste(missing, collapse = ", ")))
    }

    # validate indextype
    if(!design$indextype %in% c("FD", "FI")) {
      stop(paste("Survey design", index_idx, ": indextype must be 'FD' or 'FI'"))
    }


    # validation of n areas
    if(any(design$areas > total_areas)) {
      stop(paste("Survey design", index_idx, ": areas exceed total areas in simulation"))
    }

    # validation of parameter bounds (each should be vector of length 2)
    param_bounds <- c("q_hist_bounds", "q_proj_bounds",
                      "hyperstability_hist_bounds", "hyperstability_proj_bounds",
                      "obsError_CV_hist_bounds", "obsError_CV_proj_bounds")

    for(param in param_bounds) {
      if(length(design[[param]]) != 2 || design[[param]][2] < design[[param]][1]) {
        stop(paste("Survey design", index_idx, ":", param, "must be vector [min, max] with max >= min"))
      }
    }

    # validate survey_timing for FI surveys
    if(design$indextype == "FI") {
      if(!is.numeric(design$survey_timing) || length(design$survey_timing) != 1 ||
         design$survey_timing < 0 || design$survey_timing > 1) {
        stop(paste("Survey design", index_idx, ": survey_timing must be a single value between 0 and 1"))
      }
    }

    # for FI surveys: validation of selectivity
    # validates that selectivity indices do not exceed available selectivity objects
    if(design$indextype == "FI") {
      if(design$selectivity_hist_idx > length(IndexObj@selectivity_hist_list)) {
        stop(paste("Survey design", index_idx, ": selectivity_hist_idx exceeds available selectivity objects"))
      }
      if(design$selectivity_proj_idx > length(IndexObj@selectivity_proj_list)) {
        stop(paste("Survey design", index_idx, ": selectivity_proj_idx exceeds available selectivity objects"))
      }
    }
  } # close loop for specific surveys/CPUE validations


  #checking requirements based on index type and data required
  #check if any survey is FI or any survey uses numbers (needs N arrays)

  # index = FI or If use weight is F , stop I need N
  # for each element d in the survey_design list (any= at least one element)
  # cehck if the survey index type is FI
  # or if weigth is set to F
  # Set needs_N to TRUE if:
  # Any survey is of type "FI" (fishery-independent), or useWeight is FALSE
  needs_N <- any(sapply(IndexObj@survey_design, function(d) d$indextype == "FI")) || !IndexObj@useWeight
  if(needs_N) {
    if(is.null(simulation_result$dynamics$N)) {
      stop("Survey indices (FI) and numbers indices require N arrays. Run simulation with doDiagnostic=TRUE")
    } else {
      cat("Note: Using N data from iteration 1 for all iterations (testing mode)\n")
    }
  }

  # get the life history for calcs (need it for selectivity wrapper)
  lh <- LHwrapper(simulation_result$LifeHistoryObj, simulation_result$TimeAreaObj)

  # Results list (initialize empty list) to store multiple indices
  results_list <- list()


  # add the main loop through survey / CPUE individually

  for(index_idx in 1:n_indices) {
    design <- IndexObj@survey_design[[index_idx]]  #extract current index design

    #just to indicate what index is processing
    cat(paste("Processing", ifelse(design$indextype == "FI", "Survey", "CPUE"), index_idx,
              "- Areas:", paste(design$areas, collapse = ","), "\n"))


    # create selectivity based on the data type "FD" or "FI"
    # this section determines which  selectivity to use for calculating indices
    # in different time periods

    # historical period
    # creates historical period selectivity using fishSimGTG's historical fishery object
    if(design$indextype == "FD") {
      # for CPUE: use fishery selectivity
      # why? because that is fishsimGTG configuration
      index_selectivity_hist <- selWrapper(lh, simulation_result$TimeAreaObj,
                                           FisheryObj = simulation_result$HistFisheryObj,
                                           doPlot = FALSE)

      # For projection period: create area-specific fishery selectivities
      # initializes list for area-specific projection selectivities
      index_selectivity_proj_list <- list()

      # checks if projection fisheries exist
      if(!is.null(simulation_result$ProFisheryObj_list) &&
         length(simulation_result$ProFisheryObj_list) > 0) {
        # create selectivity for each area that this index might cover (following fishSimGTG structure)
        for(area in 1:total_areas) {
          if(area <= length(simulation_result$ProFisheryObj_list)) {
            index_selectivity_proj_list[[area]] <- selWrapper(lh, simulation_result$TimeAreaObj,
                                                              FisheryObj = simulation_result$ProFisheryObj_list[[area]],
                                                              doPlot = FALSE)
          } else {
            # uses historical if projection (ProFisheryObj_list) not available for this specific area
            index_selectivity_proj_list[[area]] <- index_selectivity_hist
          }
        }
      } else {
        # No projection fisheries exists anywhere, use historical sel for all areas
        # for example if ProFisheryObj_list is NULL or empty
        for(area in 1:total_areas) {
          index_selectivity_proj_list[[area]] <- index_selectivity_hist
        }
      }

    } else {

      # FI: use survey selectivity
      # extracts selectivity objects from provided lists and creates selectivity using sel wrappers
      hist_selectivity_obj <- IndexObj@selectivity_hist_list[[design$selectivity_hist_idx]]
      proj_selectivity_obj <- IndexObj@selectivity_proj_list[[design$selectivity_proj_idx]]

      index_selectivity_hist <- selWrapper(lh, simulation_result$TimeAreaObj,
                                           FisheryObj = hist_selectivity_obj,
                                           doPlot = FALSE)

      index_selectivity_proj <- selWrapper(lh, simulation_result$TimeAreaObj,
                                           FisheryObj = proj_selectivity_obj,
                                           doPlot = FALSE)
    }


    # initialize results for the survey/cpue
    # creates empty matrix [years × iterations] to store index values


    # similar idea to what I did with LC obs
    area_matrix <- array(NA, dim = c(years, iterations))
    rownames(area_matrix) <- paste0("Year_", 1:years)
    colnames(area_matrix) <- paste0("Iter_", 1:iterations)


    # name based on areas covered  y survey/CPUE
    if(length(design$areas) == 1) {
      matrix_name <- paste0("Area_", design$areas[1])
    } else {
      matrix_name <- paste0("Areas_", paste(design$areas, collapse = "_"))
    }
    # stores matrix in named list
    indices_by_area <- list()
    indices_by_area[[matrix_name]] <- area_matrix


    #loop through each iterations/ years combination
    for(iter in 1:iterations) {
      for(year in 1:years) {

        #check if CPUE data exists in this year
        #only calculates CPUE/Surveys for years where data is collected
        if(year %in% design$indexYears) {

          # Create a branch to account for index covering (sampling) only one area
          # or covering multiple areas

          # same patern as LC
          # checks if this survey/CPUE covers one area or multiple areas

          if(length(design$areas) == 1) {
            area <- design$areas[1] # if single area (e.g., c(1)): extract value for that area

            # get the precalc data for this single area
            if(IndexObj@useWeight && design$indextype == "FD") {
              # For CPUE biomass: get VB from this area
              area_value <- simulation_result$dynamics$VB[year, iter, area] #directly extracts vulnerable biomass from simulation results
            } else {

              # FI
              # For FI or CPUE numbers: calculate from N arrays for this area
              # manually calculates by applying selectivity to N arrays
              # and summing across GTGs.
              area_value <- 0           # starts with zero and accumulate
              # get selectivity
              # chose selectivity based on:
              # time period: historical vs projection
              # data type: FD (fishery) vs FI (survey)
              # area: For FD in proj, each area might have different fishing sel

              if(year <= historical_end) {
                selectivity <- index_selectivity_hist
              } else {
                if(design$indextype == "FD") {
                  # use area-specific fishery selectivity for projection period
                  selectivity <- index_selectivity_proj_list[[area]]  #to follow the actual configuration of the package
                } else {
                  # For FI, same selectivity for all areas this survey covers
                  selectivity <- index_selectivity_proj
                }
              }

              # sum across GTGs for this single area
              for(gtg in 1:lh$gtg) {
                N_gtg_area <- simulation_result$dynamics$N[[gtg]][, year, area]
                if(IndexObj@useWeight) {
                  survey_calc <- sum(N_gtg_area * selectivity$vul[[gtg]] * lh$W[[gtg]])
                } else {
                  survey_calc <- sum(N_gtg_area * selectivity$vul[[gtg]])
                }
                area_value <- area_value + survey_calc
              }
            }

          } else {
            # Here multiarea starts
            # For multi-area data collection: Sum data across multiple areas
            # Same approach used in length composition obs models: rowSums(true_length_array[, year, design$areas])

            if(IndexObj@useWeight && design$indextype == "FD") {
              # Sum VB across all specified areas - like rowSums() in length composition
              area_value <- sum(simulation_result$dynamics$VB[year, iter, design$areas])

            } else {
              # sum calculated values across all specified areas
              area_value <- 0

              # loop through each area and sum their contributions
              for(area in design$areas) {
                # get area-specific selectivity for this time period
                if(year <= historical_end) {
                  selectivity <- index_selectivity_hist
                } else {
                  if(design$indextype == "FD") {
                    # use area-specific fishery selectivity for projection period
                    selectivity <- index_selectivity_proj_list[[area]]
                  } else {
                    # for FI, same selectivity for all areas this survey covers
                    selectivity <- index_selectivity_proj
                  }
                }

                # calculate contribution from this area
                area_contribution <- 0
                for(gtg in 1:lh$gtg) {
                  N_gtg_area <- simulation_result$dynamics$N[[gtg]][, year, area]
                  if(IndexObj@useWeight) {
                    survey_calc <- sum(N_gtg_area * selectivity$vul[[gtg]] * lh$W[[gtg]])
                  } else {
                    survey_calc <- sum(N_gtg_area * selectivity$vul[[gtg]])
                  }
                  area_contribution <- area_contribution + survey_calc
                }
                # add this area contribution to the total (like i did in lc)
                area_value <- area_value + area_contribution
              }
            }
          }

          # get index specific parameters with bounds
          # Whether the index covers 1 area or multiple areas, it has ONE set of parameters

          # Sample parameters based on historical vs projection period
          if(year <= historical_end) {


            q_area_year <- runif(1,
                                 min = design$q_hist_bounds[1],
                                 max = design$q_hist_bounds[2])


            hyperstability_area <- runif(1,
                                         min = design$hyperstability_hist_bounds[1],
                                         max = design$hyperstability_hist_bounds[2])

            obs_CV <- runif(1,
                            min = design$obsError_CV_hist_bounds[1],
                            max = design$obsError_CV_hist_bounds[2])



          } else {
            # Projection
            q_area_year <- runif(1,
                                 min = design$q_proj_bounds[1],
                                 max = design$q_proj_bounds[2])

            hyperstability_area <- runif(1,
                                         min = design$hyperstability_proj_bounds[1],
                                         max = design$hyperstability_proj_bounds[2])

            obs_CV <- runif(1,
                            min = design$obsError_CV_proj_bounds[1],
                            max = design$obsError_CV_proj_bounds[2])
          }

          #apply the observation model to the area_value (which is either single area value or sum of multiple areas)

          if(hyperstability_area != 1) {
            index_value <- q_area_year * (area_value^hyperstability_area)
          } else {
            index_value <- q_area_year * area_value
          }

          #add observation error with bias correction
          if(obs_CV > 0) {
            obs_error <- exp(rnorm(1, mean = 0, sd = obs_CV) - 0.5 * obs_CV^2)
            index_value <- index_value * obs_error
          }

          #store the final index value in the matrix
          indices_by_area[[matrix_name]][year, iter] <- index_value
        }
      }
    }
    # store results for this survey/CPUE
    results_list[[index_idx]] <- list(
      areas = design$areas,
      indextype = design$indextype,
      indexYears = design$indexYears,
      survey_timing = if(design$indextype == "FI") design$survey_timing else NA,
      selectivity_hist_idx = if(design$indextype == "FI") design$selectivity_hist_idx else NA,
      selectivity_proj_idx = if(design$indextype == "FI") design$selectivity_proj_idx else NA,
      indices = indices_by_area
    )
  }

  # name results (name index)
  names(results_list) <- paste0(ifelse(sapply(IndexObj@survey_design, function(d) d$indextype) == "FI", "Survey_", "CPUE_"), 1:n_indices)

  # Determine the actual index type (FD, FI, or Mixed)
  actual_types <- unique(sapply(IndexObj@survey_design, function(d) d$indextype))

  if(length(actual_types) == 1) {
    # all indices are the same type
    final_indextype <- actual_types[1]  # "FD" or "FI"
  } else {
    # mixed types
    final_indextype <- "Mixed"
  }

  # return results
  return(list(
    indexID = IndexObj@indexID,
    title = IndexObj@title,
    indextype = final_indextype, # "Mixed",  # indicate this can contain both FD and FI
    useWeight = IndexObj@useWeight,
    indices = results_list,
    years_total = years,
    iterations_total = iterations,
    areas_total = total_areas,
    historical_end = historical_end,
    n_indices = n_indices
  ))
}




#---------------------------------------------------------------#
#---Step 4: simple plotting function to explore the sim index---#
#---------------------------------------------------------------#

plotIndex_simple <- function(index_result, save_plot = FALSE,
                             filename = "index_simple.jpeg") {

  years <- 1:index_result$years_total
  historical_end <- index_result$historical_end
  title <- index_result$title

  all_data <- data.frame()

  for(index_name in names(index_result$indices)) {
    index_data <- index_result$indices[[index_name]]
    indexYears <- index_data$indexYears

    for(area_name in names(index_data$indices)) {
      area_matrix <- index_data$indices[[area_name]]

      if(!is.null(area_matrix)) {
        # Median
        area_median <- apply(area_matrix, 1, median, na.rm = TRUE)
        area_median[!(years %in% indexYears)] <- NA

        median_df <- data.frame(
          year = years,
          value = area_median,
          panel = paste(index_name, area_name),
          type = "Median",
          iteration = "Median",
          stringsAsFactors = FALSE
        )
        all_data <- rbind(all_data, median_df)

        # add all iterations
        for(iter in 1:ncol(area_matrix)) {
          iter_values <- area_matrix[, iter]
          iter_values[!(years %in% indexYears)] <- NA

          iter_df <- data.frame(
            year = years,
            value = iter_values,
            panel = paste(index_name, area_name),
            type = "Iterations",
            iteration = paste0("Iter_", iter),
            stringsAsFactors = FALSE
          )
          all_data <- rbind(all_data, iter_df)
        }
      }
    }
  }

  p <- ggplot(all_data, aes(x = year, y = value)) +
    geom_line(data = subset(all_data, type == "Iterations"),
              aes(group = iteration), color = "steelblue", alpha = 0.6, size = 0.5) +

    geom_point(data = subset(all_data, type == "Iterations" & !is.na(value)),
               color = "steelblue", alpha = 0.7, size = 1) +

    geom_line(data = subset(all_data, type == "Median"),
              color = "black", size = 1.2) +
    geom_point(data = subset(all_data, type == "Median"),
               color = "black", size = 1.5) +
    facet_wrap(~ panel, scales = "free_y") +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    labs(title = title, x = "Year", y = "Index Value") +
    theme_minimal()

  if(save_plot) {
    ggsave(filename, plot = p, width = 12, height = 8, dpi = 300)
  }

  return(p)
}



# modifications
calculate_single_Index  <- function(simulation_result, IndexObj){

  # define the dimensions (to get the structure of the simulation)
  years <- dim(VB)[1]       #total years hist+future
  iterations <- dim(VB)[2]  #total iterations
  total_areas <- dim(VB)[3] #total areas
  historicalYears <- TimeAreaObj@historicalYears #total historical years- before management
  historical_end <- 1 + historicalYears

  # counts number of individual survey/CPUE programs defined in the design list.
  n_indices <- length(IndexObj@survey_design)

  if(n_indices == 0) {
    stop("survey_design list must contain at least one survey or CPUE design")
  }

  #  add validation for each survey/CPUE design
  #  lopp through each individual survey/CPUE design for validation.
  for(index_idx in 1:n_indices) {
    design <- IndexObj@survey_design[[index_idx]] #extracts current survey design from the list (for each index)

    # common required elements for surveys/CPUEs
    required_elements <- c("indextype", "areas", "indexYears",
                           "q_hist_bounds", "q_proj_bounds",
                           "hyperstability_hist_bounds", "hyperstability_proj_bounds",
                           "obsError_CV_hist_bounds", "obsError_CV_proj_bounds")

    # for FI surveys only, selectivity indices and survey_timing are needed
    if("indextype" %in% names(design) && design$indextype == "FI") {
      required_elements <- c(required_elements, "selectivity_hist_idx",
                             "selectivity_proj_idx", "survey_timing")
    }

    # check the required elements are available
    if(!all(required_elements %in% names(design))) { # check if required element exist
      missing <- required_elements[!required_elements %in% names(design)]   # identify what element is misssing
      stop(paste("Survey design", index_idx, "missing elements:", paste(missing, collapse = ", ")))
    }

    # validate indextype
    if(!design$indextype %in% c("FD", "FI")) {
      stop(paste("Survey design", index_idx, ": indextype must be 'FD' or 'FI'"))
    }


    # validation of n areas
    if(any(design$areas > total_areas)) {
      stop(paste("Survey design", index_idx, ": areas exceed total areas in simulation"))
    }

    # validation of parameter bounds (each should be vector of length 2)
    param_bounds <- c("q_hist_bounds", "q_proj_bounds",
                      "hyperstability_hist_bounds", "hyperstability_proj_bounds",
                      "obsError_CV_hist_bounds", "obsError_CV_proj_bounds")

    for(param in param_bounds) {
      if(length(design[[param]]) != 2 || design[[param]][2] < design[[param]][1]) {
        stop(paste("Survey design", index_idx, ":", param, "must be vector [min, max] with max >= min"))
      }
    }

    # validate survey_timing for FI surveys
    if(design$indextype == "FI") {
      if(!is.numeric(design$survey_timing) || length(design$survey_timing) != 1 ||
         design$survey_timing < 0 || design$survey_timing > 1) {
        stop(paste("Survey design", index_idx, ": survey_timing must be a single value between 0 and 1"))
      }
    }

    # for FI surveys: validation of selectivity
    # validates that selectivity indices do not exceed available selectivity objects
    if(design$indextype == "FI") {
      if(design$selectivity_hist_idx > length(IndexObj@selectivity_hist_list)) {
        stop(paste("Survey design", index_idx, ": selectivity_hist_idx exceeds available selectivity objects"))
      }
      if(design$selectivity_proj_idx > length(IndexObj@selectivity_proj_list)) {
        stop(paste("Survey design", index_idx, ": selectivity_proj_idx exceeds available selectivity objects"))
      }
    }
  } # close loop for specific surveys/CPUE validations


  #checking requirements based on index type and data required
  #check if any survey is FI or any survey uses numbers (needs N arrays)

  # index = FI or If use weight is F , stop I need N
  # for each element d in the survey_design list (any= at least one element)
  # cehck if the survey index type is FI
  # or if weigth is set to F
  # Set needs_N to TRUE if:
  # Any survey is of type "FI" (fishery-independent), or useWeight is FALSE
  # needs_N <- any(sapply(IndexObj@survey_design, function(d) d$indextype == "FI")) || !IndexObj@useWeight
  # if(needs_N) {
  #   if(is.null(simulation_result$dynamics$N)) {
  #     stop("Survey indices (FI) and numbers indices require N arrays. Run simulation with doDiagnostic=TRUE")
  #   } else {
  #     cat("Note: Using N data from iteration 1 for all iterations (testing mode)\n")
  #   }
  # }


  # NEW
  # Results list (initialize empty list) to store multiple indices

  #? results_list <- list()

  # Now removing the empty list and add the tibble that will be
  # filled every time step and moving up the actual_types
  # before creating the tibble:

  # Determine the actual index type (FD, FI, or Mixed)
  # then adding the calculated column  final_indextype to the tibble
  # I moved this calculation from below
  actual_types <- unique(sapply(IndexObj@survey_design, function(d) d$indextype))

  if(length(actual_types) == 1) {
    # all indices are the same type
    final_indextype <- actual_types[1]  # "FD" or "FI"
  } else {
    # mixed types
    final_indextype <- "Mixed"
  }

  # Creating the tibble
  # save the tibble in observation_return
  # each element are the columns of the tibble - creating one row at a time
  observation_return <- tibble::tibble(
    k = k,                              # current iter
    j = j,                              # current year
    indexID = IndexObj@indexID,         # index identifier
    title = IndexObj@title,             # index title
    final_indextype = final_indextype,  # "FD", "FI", or "Mixed"
    useWeight = IndexObj@useWeight,     # TRUE=biomass, FALSE=numbers
    n_indices = n_indices,              # number of indices in this index object
    years_total = years,                # total years in simulation
    iterations_total = iterations,      # total iterations in simulation
    areas_total = total_areas,          # total areas in simulation
    historical_end = historical_end     # end of historical period
  )




  # add the main loop through survey / CPUE individually

  for(index_idx in 1:n_indices) {
    design <- IndexObj@survey_design[[index_idx]]  #extract current index design

    #just to indicate what index is processing
    cat(paste("Processing", ifelse(design$indextype == "FI", "Survey", "CPUE"), index_idx,
              "- Areas:", paste(design$areas, collapse = ","), "\n"))


    # create selectivity based on the data type "FD" or "FI"
    # this section determines which  selectivity to use for calculating indices
    # in different time periods

    # historical period
    # creates historical period selectivity using fishSimGTG's historical fishery object
    if(design$indextype == "FD") {
      # for CPUE: use fishery selectivity
      # why? because that is fishsimGTG configuration
      index_selectivity_hist <- selHist[[1]]$keep # sel hist for area 1 - it is repeated for each area

      # For projection period: create area-specific fishery selectivities
      # initializes list for area-specific projection selectivities
      index_selectivity_proj_list <- selPro

    } else {

      # FI: use survey selectivity
      # extracts selectivity objects from provided lists and creates selectivity using sel wrappers
      hist_selectivity_obj <- IndexObj@selectivity_hist_list[[design$selectivity_hist_idx]]
      proj_selectivity_obj <- IndexObj@selectivity_proj_list[[design$selectivity_proj_idx]]

      index_selectivity_hist <- selWrapper(lh, TimeAreaObj,
                                           FisheryObj = hist_selectivity_obj,
                                           doPlot = FALSE)$vul

      index_selectivity_proj <- selWrapper(lh, TimeAreaObj,
                                           FisheryObj = proj_selectivity_obj,
                                           doPlot = FALSE)
    }


    # initialize results for the survey/cpue
    # creates empty matrix [years × iterations] to store index values

    #NEW:
    #remove the area_matrix

    # similar idea to what I did with LC obs
    # area_matrix <- array(NA, dim = c(years, iterations))
    # rownames(area_matrix) <- paste0("Year_", 1:years)
    # colnames(area_matrix) <- paste0("Iter_", 1:iterations)


    # # name based on areas covered  y survey/CPUE
    # if(length(design$areas) == 1) {
    #   matrix_name <- paste0("Area_", design$areas[1])
    # } else {
    #   matrix_name <- paste0("Areas_", paste(design$areas, collapse = "_"))
    # }
    # # stores matrix in named list
    # indices_by_area <- list()
    # indices_by_area[[matrix_name]] <- area_matrix


    #loop through each iterations/ years combination
    for(iter in 1:iterations) {
      for(year in 1:years) {

        #check if CPUE data exists in this year
        #only calculates CPUE/Surveys for years where data is collected
        if(j %in% design$indexYears) {

          # Create a branch to account for index covering (sampling) only one area
          # or covering multiple areas

          # same patern as LC
          # checks if this survey/CPUE covers one area or multiple areas

          if(length(design$areas) == 1) {
            area <- design$areas[1] # if single area (e.g., c(1)): extract value for that area

            # get the precalc data for this single area
            if(IndexObj@useWeight && design$indextype == "FD") {
              # For CPUE biomass: get VB from this area
              area_value <- RB[j, k, area] #directly extracts vulnerable biomass from simulation results
            } else {

              # FI
              # For FI or CPUE numbers: calculate from N arrays for this area
              # manually calculates by applying selectivity to N arrays
              # and summing across GTGs.
              area_value <- 0           # starts with zero and accumulate
              # get selectivity
              # chose selectivity based on:
              # time period: historical vs projection
              # data type: FD (fishery) vs FI (survey)
              # area: For FD in proj, each area might have different fishing sel

              if(j <= historical_end) {
                selectivity <- index_selectivity_hist
              } else {
                if(design$indextype == "FD") {
                  # use area-specific fishery selectivity for projection period
                  selectivity <- index_selectivity_proj_list[[area]]$keep  #to follow the actual configuration of the package
                } else {
                  # For FI, same selectivity for all areas this survey covers
                  selectivity <- index_selectivity_proj$vul
                }
              }

              # sum across GTGs for this single area (retained numbers)
              for(gtg in 1:lh$gtg) {
                N_gtg_area <- N[[gtg]][, j, area]
                if(IndexObj@useWeight) {
                  survey_calc <- sum(N_gtg_area * selectivity[[gtg]] * lh$W[[gtg]])
                } else {
                  survey_calc <- sum(N_gtg_area * selectivity[[gtg]])
                }
                area_value <- area_value + survey_calc
              }
            }

          } else {
            # Here multiarea starts
            # For multi-area data collection: Sum data across multiple areas
            # Same approach used in length composition obs models: rowSums(true_length_array[, year, design$areas])

            if(IndexObj@useWeight && design$indextype == "FD") {
              # Sum VB across all specified areas - like rowSums() in length composition
              area_value <- sum(RB[j, k, design$areas])

            } else {
              # sum calculated values across all specified areas
              area_value <- 0

              # loop through each area and sum their contributions
              for(area in design$areas) {
                # get area-specific selectivity for this time period
                if(j <= historical_end) {
                  selectivity <- index_selectivity_hist
                } else {
                  if(design$indextype == "FD") {
                    # use area-specific fishery selectivity for projection period
                    selectivity <- index_selectivity_proj_list[[area]]$keep
                  } else {
                    # for FI, same selectivity for all areas this survey covers
                    selectivity <- index_selectivity_proj$vul
                  }
                }

                # calculate contribution from this area
                area_contribution <- 0
                for(gtg in 1:lh$gtg) {
                  N_gtg_area <- N[[gtg]][, j, area]
                  if(IndexObj@useWeight) {
                    survey_calc <- sum(N_gtg_area * selectivity[[gtg]] * lh$W[[gtg]])
                  } else {
                    survey_calc <- sum(N_gtg_area * selectivity[[gtg]])
                  }
                  area_contribution <- area_contribution + survey_calc
                }
                # add this area contribution to the total (like i did in lc)
                area_value <- area_value + area_contribution
              }
            }
          }

          # get index specific parameters with bounds
          # Whether the index covers 1 area or multiple areas, it has ONE set of parameters

          # Sample parameters based on historical vs projection period
          if(j <= historical_end) {


            q_area_year <- runif(1,
                                 min = design$q_hist_bounds[1],
                                 max = design$q_hist_bounds[2])


            hyperstability_area <- runif(1,
                                         min = design$hyperstability_hist_bounds[1],
                                         max = design$hyperstability_hist_bounds[2])

            obs_CV <- runif(1,
                            min = design$obsError_CV_hist_bounds[1],
                            max = design$obsError_CV_hist_bounds[2])



          } else {
            # Projection
            q_area_year <- runif(1,
                                 min = design$q_proj_bounds[1],
                                 max = design$q_proj_bounds[2])

            hyperstability_area <- runif(1,
                                         min = design$hyperstability_proj_bounds[1],
                                         max = design$hyperstability_proj_bounds[2])

            obs_CV <- runif(1,
                            min = design$obsError_CV_proj_bounds[1],
                            max = design$obsError_CV_proj_bounds[2])
          }

          #apply the observation model to the area_value (which is either single area value or sum of multiple areas)

          if(hyperstability_area != 1) {
            index_value <- q_area_year * (area_value^hyperstability_area)
          } else {
            index_value <- q_area_year * area_value
          }

          #add observation error with bias correction
          if(obs_CV > 0) {
            obs_error <- exp(rnorm(1, mean = 0, sd = obs_CV) - 0.5 * obs_CV^2)
            index_value <- index_value * obs_error
          }

          #NEW

          # remove the the matrix to store the index and repleace by the tibble formatting

          # #store the final index value in the matrix
          # indices_by_area[[matrix_name]][year, iter] <- index_value

          # Adding more columns to the tibble
          # add to the tibble with appropriate column name (create column name like "CPUE_1" or "Survey_2")
          index_name <- paste0(ifelse(design$indextype == "FI", "Survey_", "CPUE_"), index_idx)
          observation_return[[index_name]] <- index_value  # add the new column to the tibble

          # Add individual index columns (from the original results_list)
          observation_return[[paste0(index_name, "_indextype")]] <- design$indextype
          observation_return[[paste0(index_name, "_areas")]] <- paste(design$areas, collapse = "_")
          observation_return[[paste0(index_name, "_indexYears")]] <- paste(design$indexYears, collapse = "_")
          observation_return[[paste0(index_name, "_survey_timing")]] <- if(design$indextype == "FI") design$survey_timing else NA
          observation_return[[paste0(index_name, "_selectivity_hist_idx")]] <- if(design$indextype == "FI") design$selectivity_hist_idx else NA
          observation_return[[paste0(index_name, "_selectivity_proj_idx")]] <- if(design$indextype == "FI") design$selectivity_proj_idx else NA
        } else {

          # For the years with no obervation (set to NA but still include) (see if(j %in% design$indexYears))
          index_name <- paste0(ifelse(design$indextype == "FI", "Survey_", "CPUE_"), index_idx)
          observation_return[[index_name]] <- NA
          observation_return[[paste0(index_name, "_indextype")]] <- design$indextype
          observation_return[[paste0(index_name, "_areas")]] <- paste(design$areas, collapse = "_")
          observation_return[[paste0(index_name, "_indexYears")]] <- paste(design$indexYears, collapse = "_")
          observation_return[[paste0(index_name, "_survey_timing")]] <- if(design$indextype == "FI") design$survey_timing else NA
          observation_return[[paste0(index_name, "_selectivity_hist_idx")]] <- if(design$indextype == "FI") design$selectivity_hist_idx else NA
          observation_return[[paste0(index_name, "_selectivity_proj_idx")]] <- if(design$indextype == "FI") design$selectivity_proj_idx else NA
        }
      } # end year
      return(observation_return)
    }
  } # end iteration
} # close the fucntion

      #End of the fucntion

    #NEW: not needed in new reporting fromat

    # store results for this survey/CPUE
    # results_list[[index_idx]] <- list(
    #   areas = design$areas,
    #   indextype = design$indextype,
    #   indexYears = design$indexYears,
    #   survey_timing = if(design$indextype == "FI") design$survey_timing else NA,
    #   selectivity_hist_idx = if(design$indextype == "FI") design$selectivity_hist_idx else NA,
    #   selectivity_proj_idx = if(design$indextype == "FI") design$selectivity_proj_idx else NA,
    #   indices = indices_by_area
    # )
  # }
  #
  # # name results (name index)
  # names(results_list) <- paste0(ifelse(sapply(IndexObj@survey_design, function(d) d$indextype) == "FI", "Survey_", "CPUE_"), 1:n_indices)

  # Determine the actual index type (FD, FI, or Mixed)
  #NEW - Moved up
  #actual_types <- unique(sapply(IndexObj@survey_design, function(d) d$indextype))
  # if(length(actual_types) == 1) {
  #   # all indices are the same type
  #   final_indextype <- actual_types[1]  # "FD" or "FI"
  # } else {
  #   # mixed types
  #   final_indextype <- "Mixed"
  # }

  # return results
#   return(list(
#     indexID = IndexObj@indexID,
#     title = IndexObj@title,
#     indextype = final_indextype, # "Mixed",  # indicate this can contain both FD and FI
#     useWeight = IndexObj@useWeight,
#     indices = results_list,
#     years_total = years,
#     iterations_total = iterations,
#     areas_total = total_areas,
#     historical_end = historical_end,
#     n_indices = n_indices
#   ))
# }
