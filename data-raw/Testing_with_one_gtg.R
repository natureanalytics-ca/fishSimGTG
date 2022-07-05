
#---------------
#Example
#---------------
#devtools::load_all()
library(fishSimGTG)
library(here)

#-----------------------------------------
#Kala
#-----------------------------------------

LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title<-"Kala"
LifeHistoryObj@speciesName<-"Naso unicornis"
LifeHistoryObj@Linf<-48
LifeHistoryObj@K<-0.43
LifeHistoryObj@t0<- -0.14
LifeHistoryObj@L50<-35.5
LifeHistoryObj@L95delta<-5.325
LifeHistoryObj@M<-0.06
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.01788
LifeHistoryObj@LW_B<-3.035
LifeHistoryObj@Steep<-0.59
LifeHistoryObj@recSD<-0 #Run with no rec var'n to see deterministic trends

HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Test"
HistFisheryObj@vulType<-"logistic"
HistFisheryObj@vulParams<-c(25.1,0.1) #Approx. knife edge
HistFisheryObj@retType<-"full"
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 1
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 100
TimeAreaObj@historicalYears = 10
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 10, ncol = 2, byrow = FALSE)

#------------------
#Higher SSB Scenario
#-----------------
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.3, 0.6)

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(25.1,0.1)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE), CPUE = c(1,2), CPUEtype='retN')

#Batch processing - 3 management strategies
stateLower<-c(35.6, 35.6, 35.6)
stateUpper<-c(0.001, 51, 200)
stateType<-c('logistic', 'slotLimit', 'slotLimit')
fileLabel<-c("Higher_option1", "Higher_option2", "Higher_option3")
projectionLabel<-c("Min", "Slot 14 - 20 inch", "Slot large")

for(sc in 1:NROW(stateLower)){

  #Size limit - changes retention, not selectivity
  ProFisheryObj@retType<-stateType[sc]
  ProFisheryObj@retParams<-c(stateLower[sc], stateUpper[sc])

  #Bag limit
  StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE), CPUE = c(1,2), CPUEtype='retN')

  #lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj, doPlot = FALSE)
  #selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = ProFisheryObj, doPlot = TRUE)


  runProjection(LifeHistoryObj = LifeHistoryObj,
                TimeAreaObj = TimeAreaObj,
                HistFisheryObj = HistFisheryObj,
                ProFisheryObj = ProFisheryObj,
                StrategyObj = StrategyObj,
                StochasticObj = StochasticObj,
                wd = here("data-test", "Kala"),
                fileName = fileLabel[sc],
                doPlot = TRUE,
                titleStrategy = projectionLabel[sc]
  )
}



relSSBscatter(wd =  here("data-test", "Kala"),
              fileName = list(
                "Higher_option1",
                "Higher_option2",
                "Higher_option3"
              ),
              facetName = c(as.list(rep("Higher biomass scenario", 3))),
              chooseArea = 0,
              proYear = 50)


relSSBseries(wd =  here("data-test", "Kala"),
             fileName = list(
               "Higher_option1",
               "Higher_option2",
               "Higher_option3"
             ),
             facetName = c(as.list(rep("Higher biomass scenario", 3))),
             chooseArea = 0,
             percentile = c(0.025, 0.975),
             doHist = TRUE,
             dpi = 300)


X<-readProjection( wd = here("data-test", "Kala"),
                  fileName = "Higher_option2"
)

Y<-readProjection( wd = here("data-test", "Kala"),
                   fileName = "Higher_option3"
)

j<-12
m<-1
k<-1
sum(sapply(1:X$lh$gtg, FUN=function(x) sum(X$dynamics$Ftotal[j,k,m]*X$selPro$keep[[x]]/(X$dynamics$Ftotal[j,k,m]*X$selPro$removal[[x]] + X$lh$LifeHistory@M)*(1-exp(-X$dynamics$Ftotal[j,k,m]*X$selPro$removal[[x]]-X$lh$LifeHistory@M))*X$dynamics$N[[x]][,j,m])))



