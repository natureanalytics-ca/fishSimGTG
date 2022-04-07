

#---------------------------------------
#Evaluate MSE
#---------------------------------------

#Roxygen header
#'Population dynamics wrapper called by runProjection
#'
#'Contains population dynamics equations. Should not be run directly, instead called by runProjection
#'
#' @param inputObject  A list of objects passed from runProjection
#' @export

evalMSE<-function(inputObject){

  #------------------
  #Unpack dataObject
  #------------------
  TimeAreaObj <- StrategyObj <- lh <- iterations <- iter <- selHist <- Ddev <- lh <- selHist <- selPro <- Cdev <- RdevMatrix <- NULL
  for(r in 1:NROW(inputObject)) assign(names(inputObject)[r], inputObject[[r]])

  controlRuleYear<-c(FALSE, rep(FALSE,(TimeAreaObj@historicalYears)), rep(TRUE, ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)))
  years <- 1 + TimeAreaObj@historicalYears + ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
  ageClasses <- lh$ageClasses
  areas <- TimeAreaObj@areas

  #--------------
  #Arrays setup
  #--------------
  SB<-array(dim=c(years, iterations, areas))
  VB<-array(dim=c(years, iterations, areas))
  catchN<-array(dim=c(years, iterations, areas))
  catchB<-array(dim=c(years, iterations, areas))
  discN<-array(dim=c(years, iterations, areas))
  discB<-array(dim=c(years, iterations, areas))
  Ftotal<-array(dim=c(years, iterations, areas))
  SPR<-array(dim=c(years, iterations))
  relSB<-array(dim=c(years, iterations))
  recN<-array(dim=c(years, iterations))

  #-----------------------------------------------
  #Setup recording of management strategy details
  #-----------------------------------------------
  decisionData<-data.frame()
  decisionAnnual<-data.frame()
  decisionLocal<-data.frame()

  #------------------------------
  #Run simulator of k iterations
  #------------------------------

  #step through iterations k
  for(k in iter[1]:iter[2]){
    #print(k)

    #-----------------------------------------
    #Initial equilibrium - year 1
    #-----------------------------------------
    is<-solveD(lh, sel = selHist, doFit = TRUE, D_type = TimeAreaObj@historicalBioType, D_in = Ddev[k])

    #Burn-in to calibrate N by area, noting effect of movement
    Ntmp <- list()
    selGroup <- 1
    yrsTmp <- (ageClasses*4)
    for (l in 1:lh$gtg){
      Ntmp[[l]]<-array(dim=c(ageClasses, yrsTmp, areas))
      for (m in 1:areas){
        Ntmp[[l]][,1,m]<-is$N[[l]]*TimeAreaObj@recArea[m]
      }
    }
    for (j in 1: yrsTmp){
      for (l in 1:lh$gtg){
        #Cohort equations + recruitment
        if(j< yrsTmp){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)
          for(m in 1:areas){
            S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=is$Feq, S_in=selHist$removal[[l]] )
            rows<-c((m-1)*dim(Ntmp[[l]])[1]+1,m*dim(Ntmp[[l]])[1])
            cols<-c(1,(dim(Ntmp[[l]])[1]*areas))
            P[rows[1]:rows[2],cols[1]:cols[2]]<- MoveMat(Surv_in=S, Move_in=TimeAreaObj@move, area_in=m)
          }
          tmp<-matrix(as.vector(Ntmp[[l]][,j,]), nrow=(dim(Ntmp[[l]])[1]*areas), ncol=1)
          Ntmp[[l]][,(j+1),]<-matrix(P%*%tmp, nrow=dim(Ntmp[[l]])[1], ncol=areas, byrow=FALSE)
          Ntmp[[l]][1,(j+1),]<- is$Req*lh$recProb[l]*TimeAreaObj@recArea
        }
      }
    }

    #Specify iniitial conditions to start sims
    N<-list()
    for(l in 1:lh$gtg){
      N[[l]]<-array(dim=c(ageClasses, years, areas))
    }
    for (l in 1:lh$gtg){
      for (m in 1:areas){
        N[[l]][,1,m]<-Ntmp[[l]][,yrsTmp,m]
        N[[l]][,2,m]<-Ntmp[[l]][,yrsTmp,m]
      }
    }

    #Arrays
    for(m in 1:areas) SB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,1,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
    SPR[1,k]<-(sum(SB[1,k,])/is$Req)/(is$B0/lh$LifeHistory@R0)
    relSB[1,k]<-sum(SB[1,k,])/is$B0
    recN[1,k]<-is$Req
    for(m in 1:areas){
      VB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist$vul[[x]]*lh$W[[x]])))
      Ftotal[1,k,m] <- is$Feq
      catchN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[1,k,m]*selHist$keep[[x]]/(Ftotal[1,k,m]*selHist$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
      catchB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[1,k,m]*selHist$keep[[x]]/(Ftotal[1,k,m]*selHist$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
      discN[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[1,k,m]*selHist$discard[[x]]/(Ftotal[1,k,m]*selHist$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
      discB[1,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[1,k,m]*selHist$discard[[x]]/(Ftotal[1,k,m]*selHist$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[1,k,m]*selHist$removal[[x]]-lh$LifeHistory@M))*N[[x]][,1,m])))
    }

    #--------------------
    #Time dynamics
    #--------------------
    for (j in 2:years){

      #Selgroup
      if(controlRuleYear[j]) selGroup <- selPro
      if(!controlRuleYear[j]) selGroup <- selHist

      #Annual regulation decisions - phase 2
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         areas = areas,
                         ageClasses = ageClasses,
                         N=N,
                         selGroup = selGroup,
                         SB=SB,
                         SPR=SPR,
                         catchN=catchN,
                         catchB=catchB,
                         Ftotal=Ftotal,
                         decisionData=decisionData,
                         decisionAnnual=decisionAnnual,
                         decisionLocal=decisionLocal,
                         Cdev=Cdev
      ),
      inputObject
      )
      if(controlRuleYear[j]) decisionAnnual<-rbind(decisionAnnual, do.call(get(StrategyObj@projectionName), list(phase=2, dataObject)))

      #Localized F at each location - phase 3
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         areas = areas,
                         ageClasses = ageClasses,
                         selGroup = selGroup,
                         N=N,
                         SB=SB,
                         SPR=SPR,
                         catchN=catchN,
                         catchB=catchB,
                         Ftotal=Ftotal,
                         decisionData=decisionData,
                         decisionAnnual=decisionAnnual,
                         decisionLocal=decisionLocal,
                         Cdev=Cdev
                        ),
                    inputObject
      )
      if(controlRuleYear[j]) { decisionLocal<-rbind(decisionLocal, do.call(get(StrategyObj@projectionName), list(phase=3, dataObject)))
      } else { decisionLocal<-rbind(decisionLocal, do.call(fixedStrategy, list(phase=3, dataObject)))}

      #SB and recruits
      for(m in 1:areas) SB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,j,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
      Rtmp<-recruit(LifeHistoryObj = lh$LifeHistory, B0=is$B0, stock=sum(SB[j,k,]))
      SPR[j,k]<-(sum(SB[j,k,])/Rtmp)/(is$B0/lh$LifeHistory@R0)
      relSB[j,k]<-sum(SB[j,k,])/is$B0
      recN[j,k]<-Rtmp*RdevMatrix[j,k]
      for (l in 1:lh$gtg) N[[l]][1,j,]<- Rtmp*lh$recProb[l]*TimeAreaObj@recArea*RdevMatrix[j,k]

      #Arrays
      for(m in 1:areas){
        VB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup$vul[[x]]*lh$W[[x]])))
        xRow<-which(decisionLocal$year==j & decisionLocal$iteration==k & decisionLocal$area==m)
        Ftotal[j,k,m] <- decisionLocal$Flocal[xRow]
        catchN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[j,k,m]*selGroup$keep[[x]]/(Ftotal[j,k,m]*selGroup$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
        catchB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[j,k,m]*selGroup$keep[[x]]/(Ftotal[j,k,m]*selGroup$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
        discN[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(Ftotal[j,k,m]*selGroup$discard[[x]]/(Ftotal[j,k,m]*selGroup$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
        discB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum(lh$W[[x]]*Ftotal[j,k,m]*selGroup$discard[[x]]/(Ftotal[j,k,m]*selGroup$removal[[x]] + lh$LifeHistory@M)*(1-exp(-Ftotal[j,k,m]*selGroup$removal[[x]]-lh$LifeHistory@M))*N[[x]][,j,m])))
      }

      #Next year abundance, move through each gtg
      for (l in 1:lh$gtg){
        if(j<years){
          P<-matrix(nrow=ageClasses*areas, ncol=ageClasses*areas)
          for(m in 1:areas){
            S<-SurvMat(ageClasses = ageClasses, M_in=lh$LifeHistory@M, F_in=Ftotal[j,k,m], S_in=selGroup$removal[[l]])
            rows<-c((m-1)*dim(N[[l]])[1]+1,m*dim(N[[l]])[1])
            cols<-c(1,(dim(N[[l]])[1]*areas))
            P[rows[1]:rows[2],cols[1]:cols[2]]<- MoveMat(Surv_in=S, Move_in=TimeAreaObj@move, area_in=m)
          }
          Ntmp<-matrix(as.vector(N[[l]][,j,]), nrow=(dim(N[[l]])[1]*areas), ncol=1)
          N[[l]][,(j+1),]<-matrix(P%*%Ntmp, nrow=dim(N[[l]])[1], ncol=areas, byrow=FALSE)
        }
      }

      #Sampling - phase 1
      if(!is.null(StrategyObj)){
        dataObject<-c(list(j=j,
                           k=k,
                           is=is,
                           areas = areas,
                           ageClasses = ageClasses,
                           N=N,
                           selGroup = selGroup,
                           SB=SB,
                           SPR=SPR,
                           catchN=catchN,
                           catchB=catchB,
                           Ftotal=Ftotal,
                           decisionData=decisionData,
                           decisionAnnual=decisionAnnual,
                           decisionLocal=decisionLocal,
                           Cdev=Cdev
        ),
        inputObject
        )
        decisionData<-rbind(decisionData, do.call(get(StrategyObj@projectionName), list(phase=1, dataObject)))
      }
    }
  }

  #save
  dynamics<-list(SB=SB, VB=VB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN)
  HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)
  return(list(dynamics=dynamics, HCR=HCR, iter=iter))
}



#---------------------------------------
#Run the projection or MSE model
#---------------------------------------

#Roxygen header
#'Run the projection or MSE model
#'
#'Function for running projections or MSE
#'
#' @param LifeHistoryObj  A LifeHistory object. Required
#' @param TimeAreaObj A TimeArea object. Required
#' @param HistFisheryObj A Fishery object that characterizes the historical dynamics. Required as it is used in initial equilibrium and historical time dynamics (if applicable)
#' @param ProFisheryObj A Fishery object used in forward projection. Optional, only used when StrategyObj is supplied
#' @param StrategyObj A Strategy object. Optional
#' @param StochasticObj A Stochastic object. Optional
#' @param wd A working directly to save output. Required
#' @param fileName A file name for output. Required
#' @param seed A value used in base::set.seed function for producing consistent set of stochastic elements. Optional
#' @param doPlot Logical whether to produce diagnostic plots upon completing simulations. Default is FALSE (no plots)
#' @param customToCluster A character vector containing name or names of custom management strategies to export to the cluster (otherwise parallel processing will fail).
#' @param titleStrategy A title for management strategy being evaluated.
#' @importFrom grDevices dev.off png rainbow
#' @importFrom graphics mtext points
#' @importFrom snowfall sfInit sfLibrary sfLapply sfRemoveAll sfStop sfExport
#' @importFrom parallel detectCores
#' @export


runProjection<-function(LifeHistoryObj, TimeAreaObj, HistFisheryObj, ProFisheryObj = NULL, StrategyObj = NULL, StochasticObj = NULL,
                        wd, fileName, seed = 1, doPlot = FALSE, customToCluster = NULL, titleStrategy = "No name"){

  #-----------------------
  #Build inputObject
  #-----------------------
  TimeAreaObj@recArea <- TimeAreaObj@recArea / sum(TimeAreaObj@recArea) #Make sure this sums to 1
  lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
  if(!is.null(lh) & lh$LifeHistory@Steep < 0.21) lh$LifeHistory@Steep <- 0.21
  if(!is.null(lh) & lh$LifeHistory@Steep > 1) lh$LifeHistory@Steep <- 1
  selHist<-selWrapper(lh, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = FALSE)
  selPro<-selWrapper(lh, TimeAreaObj, FisheryObj = ProFisheryObj, doPlot = FALSE)

  #----------------------------
  #Build stochastic components
  #----------------------------
  set.seed(seed = seed)

  #Rec devs
  RdevMatrix<-recDev(LifeHistoryObj, TimeAreaObj, StrategyObj)$Rmult

  #Initial depletion (SSB rel)
  Ddev<-bioDev(TimeAreaObj, StochasticObj)$Ddev

  #Historical cpue
  Cdev<-cpueDev(TimeAreaObj, StochasticObj)$Cdev

  #---------------------------
  #Setup parallel processing
  #---------------------------

  #Test whether we can proceed to simulations
  if(
    is.null(lh) ||
    is.null(selHist) ||
    is.null(RdevMatrix) ||
    is.null(Ddev) ||
    is.null(Cdev) ||
    length(TimeAreaObj@iterations) == 0 ||
    TimeAreaObj@iterations < 1 ||
    isTRUE(!is.null(StrategyObj) &  is.null(selPro)) ||
    isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@move)[1] | TimeAreaObj@areas != dim(TimeAreaObj@move)[2]) ||
    isTRUE(TimeAreaObj@historicalYears > 0  && isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@historicalEffort)[2] | TimeAreaObj@historicalYears != dim(TimeAreaObj@historicalEffort)[1])) ||
    isTRUE(TimeAreaObj@historicalYears + ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0) < 1) ||
    isTRUE(class(StrategyObj) == "Strategy" &&
           tryCatch({
             get(StrategyObj@projectionName)
             FALSE
            }, error = function(c) TRUE)
          ) ||
    isTRUE(class(StrategyObj) == "Strategy" && length(StrategyObj@projectionYears) == 0) ||
    isTRUE(class(StrategyObj) == "Strategy" && length(StrategyObj@projectionYears) > 0 && StrategyObj@projectionYears < 1) ||
    isTRUE(class(StrategyObj) == "Strategy" && StrategyObj@projectionName == "projectionStrategy" && length(StrategyObj@projectionParams) != 2)

  ) {
    warning("One or more components contain incomplet or erroneous information. Cannot proceed to simulation.")
    return(NULL)
  } else {

    ptm<-proc.time()
    require(snowfall)
    require(parallel)
    iterations <- floor(TimeAreaObj@iterations)

    if(detectCores()>1 & iterations >= detectCores()) {
      print("Running on multiple cores")
      cores<-min(iterations, detectCores())
      sfInit(parallel=T, cpus=cores)
      sfLibrary(fishSimGTG)
      if(!is.null(customToCluster)) sfExport(list = returnValue(customToCluster))
      input<-list()
      inputObject<-list()
      size<-floor(iterations/cores)
      for (i in 1:cores){
        input[[i]]<-c(size*(i-1)+1, ifelse(i==cores, iterations,size*i))
        inputObject[[i]]<-list(iter=c(size*(i-1)+1, ifelse(i==cores, iterations, size*i)),
                               RdevMatrix = RdevMatrix,
                               Ddev = Ddev,
                               Cdev = Cdev,
                               lh = lh,
                               selHist = selHist,
                               selPro = selPro,
                               LifeHistoryObj = LifeHistoryObj,
                               TimeAreaObj = TimeAreaObj,
                               HistFisheryObj = HistFisheryObj,
                               ProFisheryObj = ProFisheryObj,
                               StrategyObj = StrategyObj,
                               StochasticObj = StochasticObj,
                               iterations=iterations)
      }
      mseParallel<-sfLapply(inputObject, evalMSE)
      sfRemoveAll()
      sfStop()
    } else {
      mse<-evalMSE(inputObject=list(iter=c(1, iterations),
                                    RdevMatrix=RdevMatrix,
                                    Ddev=Ddev,
                                    Cdev = Cdev,
                                    lh = lh,
                                    selHist = selHist,
                                    selPro = selPro,
                                    LifeHistoryObj = LifeHistoryObj,
                                    TimeAreaObj = TimeAreaObj,
                                    HistFisheryObj = HistFisheryObj,
                                    ProFisheryObj = ProFisheryObj,
                                    StrategyObj = StrategyObj,
                                    StochasticObj = StochasticObj,
                                    iterations=iterations
                                    )
                   )
    }

    #-------------------------------
    #Ressemble from multiple cores
    #-------------------------------

    if(detectCores()>1 & iterations>=detectCores()) {

      SB<-mseParallel[[1]]$dynamics$SB
      VB<-mseParallel[[1]]$dynamics$VB
      catchB<-mseParallel[[1]]$dynamics$catchB
      catchN<-mseParallel[[1]]$dynamics$catchN
      Ftotal<-mseParallel[[1]]$dynamics$Ftotal
      discB<-mseParallel[[1]]$dynamics$discB
      discN<-mseParallel[[1]]$dynamics$discN
      SPR<-mseParallel[[1]]$dynamics$SPR
      relSB<-mseParallel[[1]]$dynamics$relSB
      recN<-mseParallel[[1]]$dynamics$recN
      decisionAnnual<-mseParallel[[1]]$HCR$decisionAnnual
      decisionLocal<-mseParallel[[1]]$HCR$decisionLocal
      decisionData<-mseParallel[[1]]$HCR$decisionData

      for (i in 2:cores){
        for(m in 1:TimeAreaObj@areas){
          SB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$SB[,input[[i]][1]:input[[i]][2],m]
          VB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$VB[,input[[i]][1]:input[[i]][2],m]
          catchB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchB[,input[[i]][1]:input[[i]][2],m]
          catchN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$catchN[,input[[i]][1]:input[[i]][2],m]
          Ftotal[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$Ftotal[,input[[i]][1]:input[[i]][2],m]
          discB[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discB[,input[[i]][1]:input[[i]][2],m]
          discN[,input[[i]][1]:input[[i]][2],m]<-mseParallel[[i]]$dynamics$discN[,input[[i]][1]:input[[i]][2],m]
        }
        SPR[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$SPR[,input[[i]][1]:input[[i]][2]]
        relSB[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$relSB[,input[[i]][1]:input[[i]][2]]
        recN[,input[[i]][1]:input[[i]][2]]<-mseParallel[[i]]$dynamics$recN[,input[[i]][1]:input[[i]][2]]
        decisionAnnual<-rbind(decisionAnnual, mseParallel[[i]]$HCR$decisionAnnual)
        decisionLocal<-rbind(decisionLocal, mseParallel[[i]]$HCR$decisionLocal)
        decisionData<-rbind(decisionData, mseParallel[[i]]$HCR$decisionData)
      }
    } else {
      SB<-mse$dynamics$SB
      VB<-mse$dynamics$VB
      catchB<-mse$dynamics$catchB
      catchN<-mse$dynamics$catchN
      Ftotal<-mse$dynamics$Ftotal
      discB<-mse$dynamics$discB
      discN<-mse$dynamics$discN
      SPR<-mse$dynamics$SPR
      relSB<-mse$dynamics$relSB
      recN<-mse$dynamics$recN
      decisionAnnual<-mse$HCR$decisionAnnual
      decisionLocal<-mse$HCR$decisionLocal
      decisionData<-mse$HCR$decisionData
    }

    #---------------
    #Save results
    #---------------

    dynamics<-list(SB=SB, VB=VB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB, recN=recN)
    HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)
    dt<-list(titleStrategy = titleStrategy, dynamics=dynamics, HCR=HCR, iterations=iterations, lh = lh,  LifeHistoryObj=LifeHistoryObj, TimeAreaObj=TimeAreaObj, HistFisheryObj=HistFisheryObj, ProFisheryObj=ProFisheryObj,  StrategyObj= StrategyObj, StochasticObj=StochasticObj)
    saveRDS(dt, file=paste(wd, "/", fileName, ".rds", sep=""))

    #--------------------------------------------------------------------------------
    #Plot results (mostly for diagnostics, these are ugly - not publication quality)
    #--------------------------------------------------------------------------------
    if(doPlot) {

      rb<-rainbow(iterations)

      #Population level plots
      png(file=paste(wd, "/", fileName, "_SPR.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$SPR[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(0,1), col=rb[1], main = "SPR")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$SPR[,k], col=rb[k])
        }
      }
      dev.off()

      png(file=paste(wd, "/", fileName, "_SBrel.png",sep=""), width=4, height=4, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$relSB[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(0,1), col=rb[1], main = "Relative spawning biomass")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$relSB[,k], col=rb[k])
        }
      }
      dev.off()

      #Recruit over time
      png(file=paste(wd, "/", fileName, "_recN.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$recN[,1], type="l", las=1, ylab="", xlab = "Year", ylim=c(min(dt$dynamics$recN),max(dt$dynamics$recN)), col=rb[1], main = "recruits N")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$recN[,k], col=rb[k])
        }
      }
      dev.off()


      #S-R
      png(file=paste(wd, "/", fileName, "_SR.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))

      is<-solveD(lh, sel = selHist, doFit = FALSE, F_in = 0.01)
      SRcurve<-t(sapply(seq(0, is$B0, length.out = 100), FUN=function(x){
        c(x/is$B0, recruit(LifeHistoryObj=dt$LifeHistoryObj, B0=is$B0, stock=x, forceR=FALSE, Rforced=0))
      }))
      plot(dt$dynamics$relSB[,1], dt$dynamics$recN[,1], type="b", las=1, ylim=c(min(dt$dynamics$recN),max(dt$dynamics$recN)), col=rb[1], ylab="Recruits", xlab = "Stock (rel SSB)", main = "Stock-recruit")
      #text(dt$dynamics$relSB[,1], dt$dynamics$recN[,1], labels=1:NROW(dt$dynamics$recN[,1]))
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$relSB[,k], dt$dynamics$recN[,k], type="b", col=rb[k])
          #text(dt$dynamics$relSB[,k], dt$dynamics$recN[,k], labels=1:NROW(dt$dynamics$recN[,k]))
        }
      }
      lines(SRcurve[,1], SRcurve[,2], type="l", col="black", las=1, ylab="Recruits", xlab = "Stock (rel SSB)", main = "Stock-recruit")

      dev.off()

      #----------------------
      # Area specific plots
      #---------------------

      #SSB
      png(file=paste0(wd, "/", fileName, "_SB_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$SB[,,m]), max(dt$dynamics$SB[,,m])), main = "Spawning biomass")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$SB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #Catch
      png(file=paste0(wd, "/", fileName, "_catchB_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$catchB[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$catchB[,,m]), max(dt$dynamics$catchB[,,m])), main = "catch in weight")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$catchB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()

      #F
      png(file=paste0(wd, "/", fileName, "_F_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$Ftotal[,1,m], type="l", las=1, ylab="", xlab = "Year", col=rb[1], ylim=c(min(dt$dynamics$Ftotal[,,m]), max(dt$dynamics$Ftotal[,,m])), main = "Fishing mortality")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$Ftotal[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #rel change in SSB
      png(file=paste0(wd, "/", fileName, "_SBchange_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m]/dt$dynamics$SB[1,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(0,2), main = "Relative spawning biomass")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$SB[,k,m]/dt$dynamics$SB[1,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()
    }
    print("Simulation time in minutes: ")
    print((proc.time()-ptm)/60)
  }
}


#---------------------------------------
#Read in data from existing projection or MSE model
#---------------------------------------

#Roxygen header
#'Read in data from existing projection or MSE model
#'
#'Function for running projections or MSE
#'
#' @param wd A working directly where output is saved. Required
#' @param fileName A file name previously used to create output. Required
#' @export

readProjection<-function(wd, fileName){
  readRDS(file=paste(wd, "/", fileName, ".rds", sep=""))
}
