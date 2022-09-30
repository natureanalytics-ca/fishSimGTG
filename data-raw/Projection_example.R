

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
LifeHistoryObj@L95delta<-5.25
LifeHistoryObj@M<-0.16
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.0136
LifeHistoryObj@LW_B<-3.109
LifeHistoryObj@Steep<-0.6
LifeHistoryObj@isHermaph<-TRUE
LifeHistoryObj@H50<-46.2
LifeHistoryObj@H95delta<-11.8
LifeHistoryObj@recSD<-0 #Run with no rec var'n to see deterministic trends
LifeHistoryObj@recRho<-0
LifeHistoryObj@R0<-10000


#---Populate a TimeArea object
#---Contains basic inputs about time and space needed to establish simulation bounds
#---The effort matrix is set as multipliers of initial equilibrium effort
TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Example"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 30
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
HistFisheryObj@vulParams<-c(40.1,0.1) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
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

#-----------------------------
#Setup stochastic object
#-----------------------------
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.3, 0.6)


#-------------------------------------------------------
#Setup fishery characteristics for the projection period
#-------------------------------------------------------

#---In this test, we are going to introduce a size limit of 35.6 cm. Ideally, size limits should be set as changes in retention, not the selectivity of the gear
#---Thus we retain the same selectivity of the historical fishery, and simply change retention
#---You should immediately notice a problem here. As selectivity does not really occur until 40 cm, while we are setting a size limit 35.6 cm (14 inches)
ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Example"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(40.1,0.1) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
ProFisheryObj@retType<-"full"
#ProFisheryObj@retParams <- c(35.6, 0.1)
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

ProFisheryObj_list<-list(ProFisheryObj, ProFisheryObj)

#------------------------
#Setup a Strategy object
#------------------------

#---The Strategy object informs the simulation that you'd like to do a projection
#---The stratgy object is used to specify effort changes (e.g. effort reduction strategies), bag limits, and spatial closures (e.g., by setting effort to 0 in a given area)
#---Since we are not modifying any of these options, we just need to create a placeholder. We will assume effort will be constant into the foreseeable future
#---The effort matrix of a Strategy object is set as a multiplier of fishing mortality in the terminal year of the historical time period - or initial eq. if no historical period specified

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99,-99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE), CPUE = c(5,8), CPUEtype = "retN", effortImpError = c(0.1, 0.2))

#----------------
#Run projection
#----------------

#---At this stage, it should be clear that we've set a size limit (even if a useless size limit) via ProFisheryObj and left all other aspects of the fishery constant
#---We will now run the projection

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              ProFisheryObj_list = ProFisheryObj_list,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test1",
              doPlot = TRUE,
              titleStrategy = "Test1",
              seed = 10
)


#---No report-ready plots are produced. However, see Batch_Projection_example for built report-ready plotting functionality
#---To create custom plots, start by reading in all the data generated by the simulation.

X<-readProjection(wd = here(),
                  fileName = "Test1"
)


