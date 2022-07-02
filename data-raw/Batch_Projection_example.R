

#---------------
#Example
#---------------
devtools::load_all()
library(fishSimGTG)
library(here)

#----------------------------
#Create a LifeHistory object
#----------------------------

#---Populate LifeHistory object
#---Contains the life history parameters
LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title<-"Hawaiian Uhu - Parrotfish"
LifeHistoryObj@speciesName<-"Chlorurus perspicillatus"
LifeHistoryObj@Linf<-53.2
LifeHistoryObj@K<-0.225
LifeHistoryObj@t0<- -1.48
LifeHistoryObj@L50<-35
LifeHistoryObj@L95<-35*1.15
LifeHistoryObj@M<-0.16
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.0136
LifeHistoryObj@LW_B<-3.109
LifeHistoryObj@Steep<-0.6
LifeHistoryObj@isHermaph<-TRUE
LifeHistoryObj@H50<-46.2
LifeHistoryObj@H95<-58
LifeHistoryObj@recSD<-0 #Run with no rec var'n to see deterministic trends

#---Populate a TimeArea object
#---Contains basic inputs about time and space needed to establish simulation bounds
#---The historical effort matrix is set as multipliers of initial equilibrium fishing mortality
#---Note that 100 iterations have been specified...this will take a few minutes to run
TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Example"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 100
TimeAreaObj@historicalYears = 10
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 10, ncol = 2, byrow = FALSE)


#---Visualize life history. Does everything make sense?
#---Optional, create a plot of life history that is useful for reports.

#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj, doPlot = TRUE)

#To save to file (for reports?)
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj, wd = here(), imageName = "LifeHistory", dpi = 300, doPlot = TRUE)

#Note that LHwrapper returns all the details of the life history
lhOut


#-----------------------------
#Setup fishery characteristics
#-----------------------------

#---Pupulate a Fishery object
#---Contains selectivity, retention and discard characteristics
#---Not sure how to set this up? Type ?selWrapper
HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Example"
HistFisheryObj@vulType<-"logistic"
HistFisheryObj@vulParams<-c(40.1,40.2) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
HistFisheryObj@retType<-"full"
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0

#---Visualize fishery vulnerability. Does everything make sense?
#---Optional, create a plot of life history that is useful for reports.

#To simply display to the console
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE)

#To save to file (for reports?)
lhOut<-LHwrapper(LifeHistoryObj, TimeAreaObj)
selWrapper(lh = lhOut, TimeAreaObj, FisheryObj = HistFisheryObj, doPlot = TRUE, wd = here(), imageName = "Vulnerability", dpi = 300)


#-----------------------
#Setup batch projections
#-----------------------

#---Here will can setup ecological scenarios, such as low or high initial relative SSB (e.g. scenarios, states of nature)
#---We can also setup multiple management options to run against each of those scenarios
#---After running these simulations, some report-ready plotting functions can be used to aggregate all of the simulations into plots.


#---Setup higher relative SSB scenario
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.3, 0.6)
StochasticObj@historicalCPUE = c(1,2) #This is used in bag limit projection
StochasticObj@historicalCPUEType = "vulN"

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Example"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(40.1,40.2) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
ProFisheryObj@retType<-"full" #We will change this as needed below in 'for' loop
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

#Batch processing - 3 management strategies
stateLmin<-c(-99, 35.6,  35.6)
stateBag<-c(2, -99, 2)
fileLabel<-c("Higher_option1", "Higher_option2", "Higher_option3")
projectionLabel<-c("Bag 2", "Min size 14 inch", "Bag 2 & min size 14 inch")

for(sc in 1:NROW(stateLmin)){

  #Size limit - changes retention, not selectivity
  if(stateLmin[sc] == -99){
    ProFisheryObj@retType<-"full"
  } else {
    ProFisheryObj@retType<-"logistic"
    ProFisheryObj@retParams<-c(stateLmin[sc],stateLmin[sc]+0.1)
  }

  #Bag limit
  StrategyObj@projectionParams<-list(bag = c(stateBag[sc], stateBag[sc]), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

  runProjection(LifeHistoryObj = LifeHistoryObj,
                TimeAreaObj = TimeAreaObj,
                HistFisheryObj = HistFisheryObj,
                ProFisheryObj = ProFisheryObj,
                StrategyObj = StrategyObj,
                StochasticObj = StochasticObj,
                wd = here(),
                fileName = fileLabel[sc],
                doPlot = TRUE,
                titleStrategy = projectionLabel[sc]
  )
}


#---Setup lower relative SSB scenario
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.3)
StochasticObj@historicalCPUE = c(1,2) #This is used in bag limit projection
StochasticObj@historicalCPUEType = "vulN"

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Example"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(40.1,40.2) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
ProFisheryObj@retType<-"full" #We will change this as needed below in 'for' loop
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

#Batch processing - 3 management strategies
stateLmin<-c(-99, 35.6,  35.6)
stateBag<-c(2, -99,  2)
fileLabel<-c("Lower_option1", "Lower_option2", "Lower_option3")
projectionLabel<-c("Bag 2", "Min size 14 inch", "Bag 2 & min size 14 inch")

for(sc in 1:NROW(stateLmin)){

  #Size limit - changes retention, not selectivity
  if(stateLmin[sc] == -99){
    ProFisheryObj@retType<-"full"
  } else {
    ProFisheryObj@retType<-"logistic"
    ProFisheryObj@retParams<-c(stateLmin[sc],stateLmin[sc]+0.1)
  }

  #Bag limit
  StrategyObj@projectionParams<-list(bag = c(stateBag[sc], stateBag[sc]), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

  runProjection(LifeHistoryObj = LifeHistoryObj,
                TimeAreaObj = TimeAreaObj,
                HistFisheryObj = HistFisheryObj,
                ProFisheryObj = ProFisheryObj,
                StrategyObj = StrategyObj,
                StochasticObj = StochasticObj,
                wd = here(),
                fileName = fileLabel[sc],
                doPlot = TRUE,
                titleStrategy = projectionLabel[sc]
  )
}


#-----------
#Charts
#-----------

#---Charts are creating using facet_wrap in ggplot2. Thus, you can group your output according to scenarios by assigning each fileName to a corresponding facetName
#---Addition options are available, try ?relSSBscatter
relSSBscatter(wd =  here(),
              fileName = list(
                "Higher_option1",
                "Higher_option2",
                "Higher_option3",
                "Lower_option1",
                "Lower_option2",
                "Lower_option3"
              ),
              facetName = c(as.list(rep("Higher biomass scenario", 3)), as.list(rep("Lower biomass scenario", 3))),
              chooseArea = 0,
              proYear = 10)

#---Addition options are available, try ?relSSBseries
relSSBseries(wd =  here(),
             fileName = list(
               "Higher_option1",
               "Higher_option2",
               "Higher_option3",
               "Lower_option1",
               "Lower_option2",
               "Lower_option3"
             ),
             facetName = c(as.list(rep("Higher biomass scenario", 3)), as.list(rep("Lower biomass scenario", 3))),
             chooseArea = 0,
             percentile = c(0.025, 0.975),
             doHist = TRUE,
             dpi = 300)




X<-readProjection( wd = here(),
                   fileName = fileLabel[1])





plot(diff(X$dynamics$recN[,1]), type = "l")

X$dynamics$Ftotal[,1,1]

plot(diff(X$dynamics$SB[,1,1]), type = "l")
