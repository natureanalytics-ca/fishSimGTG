

#---------------------------------------
#Life history wrapper for sub-cohorts
#---------------------------------------

#Roxygen header
#'Life history wrapper for sub-cohorts
#'
#'Creates the necessary age-based vectors describing life history of sub-cohorts
#'
#' @param LifeHistoryObj  A life history object.
#' @importFrom methods slot slotNames
#' @importFrom stats dnorm plogis
#' @export
#' @examples
#'LHwrapper(LifeHistoryExample)

LHwrapper<-function(LifeHistoryObj){

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
     LifeHistoryObj@L50 >= LifeHistoryObj@L95
  ) {
    NULL
  } else {

    #----------------
    #How many gtg?
    #----------------
    gtg<-13
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
    if(length(LifeHistoryObj@Tmax) == 0 || LifeHistoryObj@Tmax < 2) {
      ages<-seq(1,ceiling(-log(0.01)/LifeHistoryObj@M),1)
    } else {
      ages<-seq(1,LifeHistoryObj@Tmax,1)
    }

    #---------------------
    #Life history vectors
    #---------------------
    t0<-ifelse(length(LifeHistoryObj@t0) == 0, 0, LifeHistoryObj@t0)
    L<-lapply(gtg_Linf, FUN=function(x) x*(1-exp(-LifeHistoryObj@K*(ages-t0))))
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
      recProb = recProb,
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
#' @param LifeHistoryObj  A life history object.
#' @param FisheryObj A stock object
#' @param doProjection calculate selectivity, discard, etc. for projection time period
#' @param doPlot Creates a basic plot to visualize outcomes. Useful for ensuring parameter selections are sensible.
#' @importFrom methods slot slotNames
#' @importFrom graphics par lines legend
#' @importFrom stats median
#' @export

selWrapper<-function(LifeHistoryObj, FisheryObj, doProjection = FALSE, doPlot = FALSE){

  lh<-LHwrapper(LifeHistoryObj)

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

  #Historical
  historical<-list()
  if(is.null(lh) ||
     class(FisheryObj) != "Fishery" ||
     !(FisheryObj@historicalVulType %in%  "logistic") ||
     !(FisheryObj@historicalRetType %in%  c("full", "logistic", "slotLimit")) ||
     length(FisheryObj@historicalRetMax) == 0 ||
     FisheryObj@historicalRetMax < 0 ||
     FisheryObj@historicalRetMax > 1 ||
     length(FisheryObj@historicalDmort) == 0 ||
     FisheryObj@historicalDmort < 0 ||
     FisheryObj@historicalDmort > 1
  ) {
    historical<-NULL
  } else {
    #Vulnerability
    if(FisheryObj@historicalVulType == "logistic") {
      historical$vul<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@historicalVulParams, maxProb = 1.0))
    }

    #Retention
    if(FisheryObj@historicalRetType == "logistic") {
      historical$ret<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@historicalRetParams, maxProb = FisheryObj@historicalRetMax))
    }

    if(FisheryObj@historicalRetType == "full") {
      historical$ret<-lapply(1:lh$gtg, FUN=function(x) fullProb(L = lh$L[[x]], maxProb = FisheryObj@historicalRetMax))
    }

    if(FisheryObj@historicalRetType == "slotLimit") {
      historical$ret<-lapply(1:lh$gtg, FUN=function(x) slotProb(L = lh$L[[x]], param = FisheryObj@historicalRetParams, maxProb = FisheryObj@historicalRetMax))
    }

    #Do Keep, Dead discards, total removals
    if(
      list(NULL) %in% historical$vul ||
      list(NULL) %in% historical$ret
    ) {
      historical<-NULL
    } else {
      #Keep
      historical$keep<-lapply(1:lh$gtg, FUN=function(x) historical$vul[[x]] * historical$ret[[x]])

      #Dead discards
      historical$discard<-lapply(1:lh$gtg, FUN=function(x) historical$vul[[x]] * (1-historical$ret[[x]]) * FisheryObj@historicalDmort)

      #Removals (catch plus dead discards)
      historical$removal<-lapply(1:lh$gtg, FUN=function(x) historical$vul[[x]] * (historical$ret[[x]] + (1-historical$ret[[x]]) * FisheryObj@historicalDmort))
    }
  }

  #Projection
  projection<-list()
  if(!doProjection ||
     is.null(lh) ||
     class(FisheryObj) != "Fishery" ||
     !(FisheryObj@projectionVulType %in%  "logistic") ||
     !(FisheryObj@projectionRetType %in%  c("full", "logistic", "slotLimit")) ||
     length(FisheryObj@projectionRetMax) == 0 ||
     FisheryObj@projectionRetMax < 0 ||
     FisheryObj@projectionRetMax > 1 ||
     length(FisheryObj@projectionDmort) == 0 ||
     FisheryObj@projectionDmort < 0 ||
     FisheryObj@projectionDmort > 1
  ) {
    projection<-NULL
  } else {
    #Vulnerability
    if(FisheryObj@projectionVulType == "logistic") {
      projection$vul<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@projectionVulParams, maxProb = 1.0))
    }

    #Retention
    if(FisheryObj@projectionRetType == "logistic") {
      projection$ret<-lapply(1:lh$gtg, FUN=function(x) logisticProb(L = lh$L[[x]], param = FisheryObj@projectionRetParams, maxProb = FisheryObj@projectionRetMax))
    }

    if(FisheryObj@projectionRetType == "full") {
      projection$ret<-lapply(1:lh$gtg, FUN=function(x) fullProb(L = lh$L[[x]], maxProb = FisheryObj@projectionRetMax))
    }

    if(FisheryObj@projectionRetType == "slotLimit") {
      projection$ret<-lapply(1:lh$gtg, FUN=function(x) slotProb(L = lh$L[[x]], param = FisheryObj@projectionRetParams, maxProb = FisheryObj@projectionRetMax))
    }

    #Do Keep, Dead discards, total removals
    if(
      list(NULL) %in% projection$vul ||
      list(NULL) %in% projection$ret
    ) {
      projection<-NULL
    } else {
      #Keep
      projection$keep<-lapply(1:lh$gtg, FUN=function(x) projection$vul[[x]] * projection$ret[[x]])

      #Dead discards
      projection$discard<-lapply(1:lh$gtg, FUN=function(x) projection$vul[[x]] * (1-projection$ret[[x]]) * FisheryObj@projectionDmort)

      #Removals (catch plus dead discards)
      projection$removal<-lapply(1:lh$gtg, FUN=function(x) projection$vul[[x]] * (projection$ret[[x]] + (1-projection$ret[[x]]) * FisheryObj@projectionDmort))
    }
  }

  #doPlot
  if(doPlot){

    if(doProjection){
      par(mfcol=c(1,2), las =1)
      x<-ceiling(median(1:lh$gtg))
      plot(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$vul)[order(unlist(lh$L))], type = "l", col = "purple", lwd =3, ylim = c(0,1), xlab = "Length", ylab = "Probability", main = "Historical")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$ret)[order(unlist(lh$L))], lwd =3, col = "blue")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$keep)[order(unlist(lh$L))], lwd =3, col = "green", lty = 3)
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$discard)[order(unlist(lh$L))], lwd =3, col = "red")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$removal)[order(unlist(lh$L))], lwd =3, col = "orange", lty = 2)
      legend("topleft", legend = c("Vulnerability", "Retention", "Keep", "Dead discards", "Removals"), fill = c("purple", "blue", "green", "red", "orange"), border = "grey", bty = "n", inset=c(0, 0.1), cex = 0.8, x.intersp = 0.3)

      plot(unlist(lh$L)[order(unlist(lh$L))], unlist(projection$vul)[order(unlist(lh$L))], type = "l", col = "purple", lwd =3, ylim = c(0,1), xlab = "Length", ylab = "Probability", main = "Projection")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(projection$ret)[order(unlist(lh$L))], lwd =3, col = "blue")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(projection$keep)[order(unlist(lh$L))], lwd =3, col = "green", lty = 3)
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(projection$discard)[order(unlist(lh$L))], lwd =3, col = "red")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(projection$removal)[order(unlist(lh$L))], lwd =3, col = "orange", lty = 2)
      legend("topleft", legend = c("Vulnerability", "Retention", "Keep", "Dead discards", "Removals"), fill = c("purple", "blue", "green", "red", "orange"), border = "grey", bty = "n", inset=c(0, 0.1), cex = 0.8, x.intersp = 0.3)

    } else {

      par(mfcol=c(1,1), las = 1)
      plot(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$vul)[order(unlist(lh$L))], type = "l", col = "purple", lwd =3, ylim = c(0,1), xlab = "Length", ylab = "Probability", main = "Historical")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$ret)[order(unlist(lh$L))], lwd =3, col = "blue")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$keep)[order(unlist(lh$L))], lwd =3, col = "green", lty = 3)
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$discard)[order(unlist(lh$L))], lwd =3, col = "red")
      lines(unlist(lh$L)[order(unlist(lh$L))], unlist(historical$removal)[order(unlist(lh$L))], lwd =3, col = "orange", lty = 2)
      legend("topleft", legend = c("Vulnerability", "Retention", "Keep", "Dead discards", "Removals"), fill = c("purple", "blue", "green", "red", "orange"), border = "grey", bty = "n", inset=c(0, 0.1), cex = 0.8, x.intersp = 0.3)

    }
  }

  return(
    list(
      historical = historical,
      projection = projection
    )
  )
}



#-----------------------------------------
#Equilibrium calculations for sub-cohorts
#-----------------------------------------

#Roxygen header
#'Equilibrium conditions for sub-cohorts
#'
#'Creates the necessary age-based vectors equilibrium abundance, biomass and catch for sub-cohorts
#' @param LifeHistoryObj  A life history object.
#' @param FisheryObj A stock object
#' @param doFit Logical. When TRUE, estimates equilibrium fishing mortality based on input D_in. Ignores F_in. Default is FALSE
#' @param F_in Equilibrium fishing mortality rate. Used to calculate equilibrium conditions of the stock. Ignored when doFit = TRUE
#' @param D_type When doFit = TRUE, specifies type of equilibrium state metric that is specified in D_in (e.g., SSB depletion or SPR). Currently not in use and only SSB depletion is supported
#' @param D_in When doFit = TRUE, specifies value of equilibrium state. Currently this must be SSB depletion (value between 0 and 1)
#' @param doPlot Equilibrium length composition
#' @importFrom methods slot slotNames
#' @import ggplot2 gridExtra
#' @importFrom stats optimize
#' @export

solveD<-function(LifeHistoryObj, FisheryObj, doFit = FALSE, F_in = NULL, D_type = NULL, D_in = NULL, doPlot = FALSE){

  #------------------------------------------
  #Create life history wrapper & selectivity
  #------------------------------------------
  lh<-LHwrapper(LifeHistoryObj)
  sel<-selWrapper(LifeHistoryObj, FisheryObj, doProjection = FALSE, doPlot = FALSE)

  if(is.null(lh) ||
     is.null(sel) ||
     length(lh$LifeHistory@Steep) == 0 ||
     lh$LifeHistory@Steep <= 0.2 ||
     lh$LifeHistory@Steep > 1
  ) {
    return(NULL)
  } else {

    #----------------------------------------
    #Fitting functions
    #----------------------------------------
    min.Depletion<-function(logFmort){
      Fmort<-exp(logFmort)
      N<-lapply(1:lh$gtg, FUN=function(x) {
        Neq<-vector()
        Neq[1]<-lh$recProb[x]
        for(i in 2:lh$maxAge) Neq[i]<-exp(-lh$LifeHistory@M - Fmort*sel$historical$removal[[x]][i-1])*Neq[i-1]
        Neq
      })
      SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*lh$mat[[x]]*lh$W[[x]])))
      SPR<- SB / Wbar
      D<-(4*lh$LifeHistory@Steep*SPR+lh$LifeHistory@Steep-1)/(5*lh$LifeHistory@Steep-1)
      (D-D_in)^2
    }

    #---------------
    #Wbar
    #---------------
    N<-lapply(1:lh$gtg, FUN=function(x) {
      Neq<-vector()
      Neq[1]<-lh$recProb[x]
      for(i in 2:lh$maxAge) Neq[i]<-exp(-lh$LifeHistory@M)*Neq[i-1]
      Neq
    })
    Wbar<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*lh$mat[[x]]*lh$W[[x]])))

    #-------------
    #Get Feq
    #-------------

    ##FIX PROBLEM WITH optimize

    if(doFit){
      Feq<-exp(optimize(min.Depletion, lh, sel, D_in, Wbar, lower=-14, upper=1.1, maximum=FALSE, tol=0.00000001)$minimum)
    } else {
      Feq<-F_in
    }

    #-------------------------------------------------------------------------------------------
    #Get current SSB depletion
    #Re-cacl regardless of having D_in because in some cases fit cannot deplete stock to target
    #-------------------------------------------------------------------------------------------
    N<-lapply(1:lh$gtg, FUN=function(x) {
      Neq<-vector()
      Neq[1]<-lh$recProb[x]
      for(i in 2:lh$maxAge) Neq[i]<-exp(-lh$LifeHistory@M - Feq*sel$historical$removal[[x]][i-1])*Neq[i-1]
      Neq
    })
    YPR<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq*sel$historical$keep[[x]]/(Feq*sel$historical$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq*sel$historical$removal[[x]]-lh$LifeHistory@M))*N[[x]])))
    SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*lh$mat[[x]]*lh$W[[x]])))
    SPR<-SB / Wbar
    D<-max(0, (4*lh$LifeHistory@Steep*SPR+lh$LifeHistory@Steep-1)/(5*lh$LifeHistory@Steep-1))

    #------------------------------
    #Scale stock size and catches
    #------------------------------
    Req<-recruit(LifeHistoryObj = lh$LifeHistory, B0=Wbar, stock=Wbar*D, forceR=FALSE)
    N<-lapply(1:lh$gtg, FUN=function(x) {
      Neq<-vector()
      Neq[1]<-lh$recProb[x]*Req
      for(i in 2:lh$maxAge) Neq[i]<-exp(-lh$LifeHistory@M - Feq*sel$historical$removal[[x]][i-1])*Neq[i-1]
      Neq
    })
    SB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*lh$mat[[x]]*lh$W[[x]])))
    VB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]]*sel$historical$vul[[x]]*lh$W[[x]])))
    catchN<-sum(sapply(1:lh$gtg, FUN=function(x) sum(Feq*sel$historical$keep[[x]]/(Feq*sel$historical$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq*sel$historical$removal[[x]]-lh$LifeHistory@M))*N[[x]])))
    catchB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq*sel$historical$keep[[x]]/(Feq*sel$historical$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq*sel$historical$removal[[x]]-lh$LifeHistory@M))*N[[x]])))
    discN<-sum(sapply(1:lh$gtg, FUN=function(x) sum(Feq*sel$historical$discard[[x]]/(Feq*sel$historical$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq*sel$historical$removal[[x]]-lh$LifeHistory@M))*N[[x]])))
    discB<-sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Feq*sel$historical$discard[[x]]/(Feq*sel$historical$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Feq*sel$historical$removal[[x]]-lh$LifeHistory@M))*N[[x]])))
    B0<-Wbar*lh$LifeHistory@R0


    if(doPlot) {



      tmp1<-data.frame(
        l = unlist(lapply(1:lh$gtg, FUN=function(x){
        rep(lh$L[[x]], N[[x]])
      })))


      tmp2<-data.frame(
        l = unlist(lapply(1:lh$gtg, FUN=function(x){
        rep(lh$L[[x]], N[[x]]*sel$historical$vul[[x]])
      })))

      p1<-ggplot(tmp1, aes(x=~l)) +
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

      p2<-ggplot(tmp2, aes(x=~l)) +
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


