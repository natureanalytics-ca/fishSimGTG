

#---------------
#Example
#---------------
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
#---The effort matrix is set as multipliers of initial equilibrium effort
TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Example"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 3
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 50, ncol = 2, byrow = FALSE)


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


#--------------------------------------------------
#At this stage, we can simulate historical dynamics
#--------------------------------------------------

#---Notice that the item TimeAreaObj@historicalEffort was set for constant effort for 10 years
#---The output is constant SSB and catch because we specified the same effort conditions to occur into the future
#---The runProjection function is the master function for conducting simulation modeling
#---The runProjection funtion produces many low quality plots. These are meant for the analyst to verify the results.
runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              wd = here(),
              fileName = "Test1",
              doPlot = TRUE,
              titleStrategy = "Test1"
)

#---What if we change the historical effort trend?
#---An increasing trend that doubles effort after 10 years
TimeAreaObj@historicalEffort<-matrix(seq(1,3,length.out = 50), nrow = 50, ncol = 2, byrow = FALSE)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              wd = here(),
              fileName = "Test2",
              doPlot = TRUE,
              titleStrategy = "Test2"
)

#---Now, lets add uncertainty in the initial depletion
#---Current, we have the following parameters specificed TimeAreaObj@historicalBio & TimeAreaObj@historicalBioType, which specify a constant SSB/SSB0 for all iterations
#---But when we use a Stochastic object, it overrides this constant value
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.3, 0.6)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test3",
              doPlot = TRUE,
              titleStrategy = "Test3"
)

#---Finally, lets introduce inter-annual recruitment variation
#---Recruitment variation is specified in the LifeHistory object.
LifeHistoryObj@recSD<-0.6
runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = HistFisheryObj,
              StochasticObj = StochasticObj,
              wd = here(),
              fileName = "Test4",
              doPlot = TRUE,
              titleStrategy = "Test4"
)
