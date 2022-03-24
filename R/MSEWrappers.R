

############################################################
#Wrappers
############################################################
evalMSE<-function(inputObject){

  #Unpack dataObject
  for(r in 1:NROW(inputObject)) assign(names(inputObject)[r], inputObject[[r]])

  #iter
  #RdevMatrix = RdevMatrix,
  #Ddev = Ddev,
  #Cdev = Cdev,
  #lh = lh,
  #selHist = selHist,
  #selPro = selPro,
  #LifeHistoryObj = LifeHistoryObj,
  #TimeAreaObj = TimeAreaObj,
  #FisheryObj = FisheryObj,
  #StrategyObj = StrategyObj,
  #iterations=iterations

  controlRuleYear<-c(rep(FALSE,(TimeAreaObj@historicalYears)), rep(TRUE, ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)))
  years <- TimeAreaObj@historicalYears + ifelse(class(StrategyObj) == "Strategy"  && length(StrategyObj@projectionYears) > 0, StrategyObj@projectionYears, 0)
  ageClasses <- lh$ageClasses
  areas <- TimeAreaObj@areas

  #ref<-refInfo(iterations, iter, lh)

  #Arrays
  SB<-array(dim=c(years, iterations, areas))
  VB<-array(dim=c(years, iterations, areas))
  catchN<-array(dim=c(years, iterations, areas))
  catchB<-array(dim=c(years, iterations, areas))
  discN<-array(dim=c(years, iterations, areas))
  discB<-array(dim=c(years, iterations, areas))
  Ftotal<-array(dim=c(years, iterations, areas))
  SPR<-array(dim=c(years, iterations))
  relSB<-array(dim=c(years, iterations))

  #Setup recording of management strategy
  decisionData<-data.frame()
  decisionAnnual<-data.frame()
  decisionLocal<-data.frame()

  ##############################
  #Run temporal dynamics
  ##############################
  #step through iterations k
  for(k in iter[1]:iter[2]){
    #print(k)

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
      }
    }

    for (j in 1:years){

      #Annual regulation decisions - phase 2
      dataObject<-list(j=j, k=k, is=is, areas = areas, ageClasses = ageClasses, N=N, SB=SB, SPR=SPR, catchN=catchN, catchB=catchB, Ftotal=Ftotal, decisionData=decisionData, decisionAnnual=decisionAnnual, decisionLocal=decisionLocal, Cdev=Cdev)

      #FIX THIS paramMSE
      if(controlRuleYear[j]) decisionAnnual<-rbind(decisionAnnual, do.call(paramMSE$mseName, list(phase=2, dataObject)))

      #Localized F at each location - phase 3
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         areas = areas,
                         ageClasses = ageClasses,
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
      #FIX THIS
      if(controlRuleYear[j]) { decisionLocal<-rbind(decisionLocal, do.call(paramMSE$mseName, list(phase=3, dataObject)))
      } else { decisionLocal<-rbind(decisionLocal, do.call(fixedRule, list(phase=3, dataObject)))}

      #Selgroup
      if(controlRuleYear[j]) selGroup <- selPro
      if(!controlRuleYear[j]) selGroup <- selHist

      #SB and recruits
      for(m in 1:areas) SB[j,k,m] <- sum(sapply(1:lh$gtg, FUN=function(x) sum((N[[x]][,j,m]*lh$mat[[x]]*lh$W[[x]])[2:ageClasses])))
      Rtmp<-recruit(LifeHistoryObj = lh$LifeHistory, B0=is$B0, stock=sum(SB[j,k,]))
      SPR[j,k]<-(sum(SB[j,k,])/Rtmp)/(is$B0/lh$LifeHistory@R0)
      relSB[j,k]<-sum(SB[j,k,])/is$B0
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
      dataObject<-c(list(j=j,
                         k=k,
                         is=is,
                         areas = areas,
                         ageClasses = ageClasses,
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
      #FIX THIS decisionData<-rbind(decisionData, do.call(paramMSE$mseName, list(phase=1, dataObject)))
    }
  }

  #save
  dynamics<-list(SB=SB, VB=VB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB)
  HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)
  return(list(dynamics=dynamics, HCR=HCR, iter=iter))
}





runProjection<-function(LifeHistoryObj, TimeAreaObj, HistFisheryObj, ProFisheryObj = NULL, StrategyObj = NULL, StochasticObj = NULL,  wd, fileName, seed = 1, doPlot = FALSE){

  #-----------------------
  #Build inputObject
  #-----------------------
  TimeAreaObj@recArea <- TimeAreaObj@recArea / sum(TimeAreaObj@recArea) #Make sure this sums to 1
  lh<-LHwrapper(LifeHistoryObj, TimeAreaObj)
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
    isTRUE(TimeAreaObj@areas != dim(TimeAreaObj@historicalEffort)[2] | TimeAreaObj@historicalYears != dim(TimeAreaObj@historicalEffort)[1])
  ) {
    warning("One or more components are incomplete or contain erroneous information. Cannot proceed to simulation.")
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
      #sfLibrary(truncnorm)
      sfLibrary(fishSimGTG)
      #sfExportAll()
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
                               FisheryObj = FisheryObj,
                               StrategyObj = StrategyObj,
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
                                    FisheryObj = FisheryObj,
                                    StrategyObj = StrategyObj,
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
      decisionAnnual<-mseParallel[[1]]$HCR$decisionAnnual
      decisionLocal<-mseParallel[[1]]$HCR$decisionLocal
      decisionData<-mseParallel[[1]]$HCR$decisionData

      for (i in 2:cores){
        for(m in 1:paramObjects$areas){
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
      decisionAnnual<-mse$HCR$decisionAnnual
      decisionLocal<-mse$HCR$decisionLocal
      decisionData<-mse$HCR$decisionData
    }

    #---------------
    #Save results
    #---------------

    dynamics<-list(SB=SB, VB=VB, catchB=catchB, catchN=catchN, Ftotal=Ftotal, discB=discB, discN=discN, SPR=SPR, relSB=relSB)
    HCR<-list(decisionLocal=decisionLocal, decisionAnnual=decisionAnnual, decisionData=decisionData)
    dt<-list(dynamics=dynamics, HCR=HCR, iterations=iterations, lh = lh,  LifeHistoryObj=LifeHistoryObj, TimeAreaObj=TimeAreaObj, HistFisheryObj=HistFisheryObj, ProFisheryObj=ProFisheryObj,  StrategyObj= StrategyObj, StochasticObj=StochasticObj)
    saveRDS(list(mse=dt), file=paste(wd, "/", fileName, "_MSE_PERFORMANCE.rda", sep=""))

    #--------------------------------------------------------------------------------
    #Plot results (mostly for diagnostics, these are ugly - not publication quality)
    #--------------------------------------------------------------------------------
    if(doPlot) {

      rb<-rainbow(iterations)

      #Population level plots
      png(file=paste(wd, "/", fileName, "_MSE_SPR.png",sep=""), width=4, height=4, units="in", res=300, bg="white", pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$SPR[,1], type="l", las=1, ylab="", ylim=c(0,1), col=rb[1], main = "SPR")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$SPR[,k], col=rb[k])
        }
      }
      dev.off()

      png(file=paste(wd, "/", fileName, "_MSE_SBrel.png",sep=""), width=4, height=4, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      plot(dt$dynamics$relSB[,1], type="l", las=1, ylab="", ylim=c(0,1), col=rb[1], main = "Relative spawning biomass")
      if(iterations > 1){
        for(k in 2:iterations){
          lines(dt$dynamics$relSB[,k], col=rb[k])
        }
      }
      dev.off()

      #----------------------
      # Area specific plots
      #---------------------

      #SSB
      png(file=paste0(wd, "/", fileName, "_MSE_SB_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(min(dt$dynamics$SB[,,m]), max(dt$dynamics$SB[,,m])), main = "Spawning biomass")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$SB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #Catch
      png(file=paste0(wd, "/", fileName, "_MSE_catchB_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$catchB[,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(min(dt$dynamics$catchB[,,m]), max(dt$dynamics$catchB[,,m])), main = "catch in weight")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$catchB[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()

      #F
      png(file=paste0(wd, "/", fileName, "_MSE_F_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$Ftotal[,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(min(dt$dynamics$Ftotal[,,m]), max(dt$dynamics$Ftotal[,,m])), main = "Fishing mortality")
        mtext(paste("Area", m), side=3, font=2, line=0.1, adj=0)
        if(iterations > 1){
          for(k in 2:iterations){
            lines(dt$dynamics$Ftotal[,k,m], type="l", las=1, ylab="",col=rb[k])
          }
        }
      }
      dev.off()


      #rel change in SSB
      png(file=paste0(wd, "/", fileName, "_MSE_SBchange_Area.png"), width=9, height=7, units="in", res=300, bg="white",pointsize=12)
      par(mfrow=c(ceiling(dt$TimeAreaObj@areas/2+0.5),2), mar=c(4,4,3,1))
      for(m in 1:dt$TimeAreaObj@areas) {
        plot(dt$dynamics$SB[,1,m]/dt$dynamics$SB[1,1,m], type="l", las=1, ylab="", col=rb[1], ylim=c(0,1.2), main = "Relative spawning biomass")
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
