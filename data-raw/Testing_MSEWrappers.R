
###################################
#Historical dynamics
###################################

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
LifeHistoryObj@Steep <- 0.6

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
StochasticObj@historicalBio = c(0.2, 0.8)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)


#------------------------------------
#Test of rebuilding to unfished level
#------------------------------------

#Turn off fishing, does stock return to B0?
devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 1.2

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(0:0, nrow=50, ncol=2)

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.9)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)



#------------------------------------
#Create a trend in historical effort
#------------------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.7

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
#TimeAreaObj@historicalEffort<-matrix(seq(1,2,length.out = 50), nrow=50, ncol=2, byrow = FALSE)
TimeAreaObj@historicalEffort<-matrix(c(rep(0.1, 50)), nrow=50, ncol=2, byrow = FALSE)

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.9)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)



#------------------------------------
#Introduce recruitment variation
#------------------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0.6
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.7

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow=50, ncol=2, byrow = FALSE)

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.9)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)


#------------------------------------
#Introduce recruitment variation + F trend
#------------------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0.6
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.7

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(seq(1,2,length.out = 50), nrow=50, ncol=2, byrow = FALSE)

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.9)

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)

#--------------------------------
#Example using prototype inputs
#-------------------------------

devtools::load_all()
library(here)

LifeHistoryObj<- new("LifeHistory")
FisheryObj<-new("Fishery")
TimeAreaObj<-new("TimeArea")
StochasticObj<-new("Stochastic")

runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, StochasticObj=StochasticObj, wd = here("data-test"), fileName = "constantModel", doPlot = TRUE)


#####################################
#Projection function testing
#####################################


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
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
#TimeAreaObj@historicalEffort<-matrix(1:1, nrow=50, ncol=2)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2))


runProjection(LifeHistoryObj = LifeHistoryObj, TimeAreaObj = TimeAreaObj, HistFisheryObj = FisheryObj, ProFisheryObj = FisheryObj, StrategyObj = StrategyObj, wd = here("data-test"), fileName = "projectionModel", doPlot = TRUE)



#-----------------------------------------------------
#Constant dynamics test w multiple starting conditions
#-------------------------------------------------

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(1:1, nrow=50, ncol=2))

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = FisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
              )

#------------------------------------
#Test of rebuilding to unfished level
#------------------------------------

#Turn off fishing, does stock return to B0?


devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.4

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(0:0, nrow=50, ncol=2))

StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = FisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)

dt<-readRDS(file = here("data-test", "projectionModel_MSE.rds"))


#---------------------------
#Create a trend in effort
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.4

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
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
#StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(c(seq(1,2,length.out = 50), rep(2, 50)), nrow=100, ncol=2, byrow = FALSE))
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(rep(0.1, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = FisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)

#---------------------------
#Create a change in size limit
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.4

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
#TimeAreaObj@historicalEffort<-matrix(1:1, nrow=50, ncol=2)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(-99, -99), effort = matrix(rep(1, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)

#---------------------------
#Create a change in bag limit
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(50,75)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(1, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)


#---------------------------
#Create a change in bag limit + size limit
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(1, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)

#---------------------------
#Create a change in bag limit + size limit + effort reduction
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 0
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix()


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(0.5, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)


#---------------------------
#Create a change in bag limit + size limit + effort reduction w/ historical period
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 5
TimeAreaObj@historicalYears = 50
TimeAreaObj@historicalBio = 0.5
TimeAreaObj@historicalBioType = "relB"
TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 50, ncol = 2, byrow = FALSE)


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(0.5, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)


#---------------------------
#Create a change in bag limit + size limit + effort reduction w/ historical period + plus rec var'n
#---------------------------

devtools::load_all()
library(here)

LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0.6
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

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
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 50, ncol = 2, byrow = FALSE)


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(0.5, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"

runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE
)


#---------------------------
#Create custom MP
#---------------------------

devtools::load_all()
library(here)
library(fishSimGTG)
LifeHistoryObj <- LifeHistoryExample
LifeHistoryObj@recSD <- 0.6
LifeHistoryObj@recRho <- 0
LifeHistoryObj@Steep <- 0.8

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@vulType<-"logistic"
FisheryObj@vulParams<-c(50,75)
FisheryObj@retType<-"full"
FisheryObj@retMax <- 1
FisheryObj@Dmort <- 0

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(75,85)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

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
TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 50, ncol = 2, byrow = FALSE)


StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 100
StrategyObj@projectionName<-"projectionStrategy_NEW"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(rep(0.5, 100), nrow=100, ncol=2, byrow = FALSE))


StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.2, 0.8)
StochasticObj@historicalCPUE = c(5,5)
StochasticObj@historicalCPUEType = "vulN"


projectionStrategy_NEW<-function(phase, dataObject){

  #Unpack dataObject
  j <- TimeAreaObj <- areas <- StrategyObj <- is <- k <- StochasticObj <- lh <- N <- selHist <- Cdev <- selGroup <- NULL
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])

  #Booking keeping for year for items in StrategyObj
  yr <- j - TimeAreaObj@historicalYears - 1

  if(phase==3){

    Flocal<-data.frame()
    for (m in 1:areas){

      bag <- StrategyObj@projectionParams[[1]][m]

      if(bag == -99){

        #Apply Flocal
        Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*is$Feq
        Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

      } else {

        if(StochasticObj@historicalCPUEType == "vulN") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist$keep[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup$keep[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*is$Feq*Fmult
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }

        if(StochasticObj@historicalCPUEType == "vulB") {

          #Initial equilibrium vulnerable N
          Nvul<-sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,1,m]*selHist$keep[[x]]*lh$W[[x]])))

          #Specify assumed initial lambda
          lambdaInitial <- Cdev[k]

          #Solove for q
          q<-lambdaInitial/Nvul

          #Get current F multiplier
          lambda<-q*sum(sapply(1:lh$gtg, FUN=function(x) sum(N[[x]][,j,m]*selGroup$keep[[x]]*lh$W[[x]])))
          probs<-c(dpois(0:(bag-1), lambda), 1-ppois(bag-1,lambda))
          probs<-probs/sum(probs)
          nm<-sum(0:bag*probs)
          Fmult<-min(nm/lambda, 1.0)

          #Apply Flocal
          Ftmp<-StrategyObj@projectionParams[[2]][yr,m]*is$Feq*Fmult
          Flocal<-rbind(Flocal, c(j, k, m, Ftmp))

        }
      }
    }
    return(list(year=Flocal[,1], iteration=Flocal[,2], area=Flocal[,3],  Flocal=Flocal[,4]))
  }
}


runProjection(LifeHistoryObj = LifeHistoryObj,
              TimeAreaObj = TimeAreaObj,
              HistFisheryObj = FisheryObj,
              ProFisheryObj = ProFisheryObj,
              StrategyObj = StrategyObj,
              StochasticObj = StochasticObj,
              wd = here("data-test"),
              fileName = "projectionModel",
              doPlot = TRUE,
              customToCluster = "projectionStrategy_NEW"
)

