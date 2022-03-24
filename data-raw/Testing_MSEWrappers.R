
#------------------------
#Constant dynamics test
#------------------------

#Do initial eq. conditions carry forward if F is constant?
#Can we achieve this with the minimum required inputs?

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.5

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 1
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow=50, ncol=2)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)


#------------------------
#Constant dynamics test 2
#------------------------

#Do initial eq. conditions carry forward if F is constant?
#Add replication at variety of depletion levels

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.5

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 10
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow=50, ncol=2)

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.9)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)






StochasticObj@historicalCPUE = c(5, 10)

HistFisheryObj<-FisheryObj


set.seed(1)
RdevMatrix<-recDev(LifeHistoryObj, TimeAreaObj, StrategyObj = NULL)$Rmult
RdevMatrix

Ddev<-bioDev(TimeAreaObj, StochasticObj = StochasticObj)$Ddev
Ddev

Cdev<-cpueDev(TimeAreaObj, StochasticObj = StochasticObj)$Cdev
Cdev
