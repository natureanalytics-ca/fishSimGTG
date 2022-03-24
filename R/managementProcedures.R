


#----------------------------------------------------------------------
#Fixed rule - for fixed harvest time period (not for MSE time period)
#----------------------------------------------------------------------

#Roxygen header
#'Historical fishing pressure
#'
#' @param phase Management procedures are coded in three phases: 1 - data collection, 2 - a decision making process, 3 - conversion of that process into annual F
#' @param dataObject The needed inputs to the management procedure
#' @export

fixedRule<-function(phase, dataObject){
  #Unpack dataObject
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==3){
    #Create a temp data frame of fishing mortalities by area
    Flocal<-data.frame()
    for (m in 1:areas) Flocal<-rbind(Flocal, c(j, k, m, TimeAreaObj@historicalEffort[j,m]*is$Feq))
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}


#-------------------------------------------------------------------
#Projection rule - no harvest control rule just projections
#-------------------------------------------------------------------

constantF_bag<-function(phase, dataObject){

  #Unpack dataObject
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  if(phase==3){

    Flocal<-data.frame()
    Dlocal<-data.frame()

    for(l in 1:paramObjects$stocks){

      for (m in 1:paramObjects$areas){

        bag<-paramMSE$bag[l,m]

        if(bag == -99){

          #Apply Flocal
          Ftmp<-is$initialMat[[l]]$Feq*paramMSE$fisher[j,m]
          Flocal<-rbind(Flocal, c(j, k, l, m, Ftmp))

          #Apply Dlocal
          Dlocal<-rbind(Dlocal, c(j, k, l, m, 0.0))

        } else {
          #Initial vulnerable N
          Nvul<-sum(N[[l]][,1,k,m]*ref$Sel[[l]][,1,k,m])

          #Specify assumed initial lamdba
          lambdaInitial <- cpue[k,l]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(N[[l]][,j,k,m]*ref$Sel[[l]][,2,k,m])

          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)

          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-is$initialMat[[l]]$Feq*paramMSE$fisher[j,m]*Fmult
          Flocal<-rbind(Flocal, c(j, k, l, m, Ftmp))

          #Apply Dlocal
          d<- -log(1-paramMSE$Dmort[l])*paramInitial$fracCatch[l]
          Dtmp<- is$initialMat[[l]]$Feq*paramMSE$fisher[j,m]*(1-Fmult)*d
          Dlocal<-rbind(Dlocal, c(j, k, l, m, Dtmp))
        }


      }
    }
    return(list(year=Flocal[,1], iteration=Flocal[,2], stock=Flocal[,3], area=Flocal[,4],  Flocal=Flocal[,5], Dlocal=Dlocal[,5]))
  }
}
