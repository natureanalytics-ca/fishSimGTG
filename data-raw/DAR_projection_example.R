


#library(fishSimGTG)
devtools::load_all()
library(here)

#-----------------------------------------
#Kole
#-----------------------------------------

LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title<-"Kole"
LifeHistoryObj@speciesName<-"Ctenochaetus strigosus"
LifeHistoryObj@Linf<-17.7
LifeHistoryObj@K<-0.423
LifeHistoryObj@t0<- -0.51
LifeHistoryObj@L50<-8.4
LifeHistoryObj@L95<-8.4*1.15
LifeHistoryObj@M<-0.08
LifeHistoryObj@L_type<-"FL"
LifeHistoryObj@L_units<-"cm"
LifeHistoryObj@LW_A<-0.046
LifeHistoryObj@LW_B<-2.85
LifeHistoryObj@Steep<-0.54
LifeHistoryObj@recSD<-0 #Run with no rec var'n to see deterministic trends

HistFisheryObj<-new("Fishery")
HistFisheryObj@title<-"Test"
HistFisheryObj@vulType<-"logistic"
HistFisheryObj@vulParams<-c(10.2,10.3) #Approx. knife edge
HistFisheryObj@retType<-"full"
HistFisheryObj@retMax <- 1
HistFisheryObj@Dmort <- 0

TimeAreaObj<-new("TimeArea")
TimeAreaObj@title = "Test"
TimeAreaObj@gtg = 13
TimeAreaObj@areas = 2
TimeAreaObj@recArea = c(0.99, 0.01)
TimeAreaObj@iterations = 10
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
StochasticObj@historicalCPUE = c(7,11)
StochasticObj@historicalCPUEType = "vulN"

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(10.2,10.3)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

#Batch processing - 3 management strategies
stateLmin<-c(10.2, 12.7,  12.7)
stateBag<-c(20, -99, 20)
fileLabel<-c("Higher_option1", "Higher_option2", "Higher_option3")
projectionLabel<-c("Bag 20", "Min size 5 inch", "Bag 20 & min size 5 inch")

for(sc in 1:NROW(stateLmin)){

  #Size limit
  ProFisheryObj@vulParams<-c(stateLmin[sc],stateLmin[sc]+0.1)

  #Bag limit
  StrategyObj@projectionParams<-list(bag = c(stateBag[sc], stateBag[sc]), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

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


#------------------
#Lower SSB Scenario
#-----------------
StochasticObj<-new("Stochastic")
StochasticObj@historicalBio = c(0.1, 0.3)
StochasticObj@historicalCPUE = c(7,11)
StochasticObj@historicalCPUEType = "vulN"

ProFisheryObj<-new("Fishery")
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(10.2,10.3)
ProFisheryObj@retType<-"full"
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0

StrategyObj <- new("Strategy")
StrategyObj@projectionYears <- 50
StrategyObj@projectionName<-"projectionStrategy"
StrategyObj@projectionParams<-list(bag = c(5, 5), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

#Batch processing - 3 management strategies
stateLmin<-c(10.2, 12.7,  12.7)
stateBag<-c(20, -99, 20)
fileLabel<-c("Lower_option1", "Lower_option2", "Lower_option3")
projectionLabel<-c("Bag 20", "Min size 5 inch", "Bag 20 & min size 5 inch")

for(sc in 1:NROW(stateLmin)){

  #Size limit
  ProFisheryObj@vulParams<-c(stateLmin[sc],stateLmin[sc]+0.1)

  #Bag limit
  StrategyObj@projectionParams<-list(bag = c(stateBag[sc], stateBag[sc]), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE))

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

#-----------
#Charts
#-----------

relSSBscatter(wd =  here("data-test", "Kala"),
                        fileName = list(
                          "Higher_option1",
                          "Higher_option2",
                          "Higher_option3",
                          "Lower_option1",
                          "Lower_option2",
                          "Lower_option3",
                          "Higher_option1",
                          "Higher_option2",
                          "Higher_option3"
                        ),
                        facetName = c(as.list(rep("Higher biomass scenario", 3)), as.list(rep("Lower biomass scenario", 3)), as.list(rep("Other biomass scenario", 3))),
                        chooseArea = 0,
                        proYear = 10)


relSSBseries(wd =  here("data-test", "Kala"),
             fileName = list(
               "Higher_option1",
               "Higher_option2",
               "Higher_option3",
               "Lower_option1",
               "Lower_option2",
               "Lower_option3",
               "Higher_option1",
               "Higher_option2",
               "Higher_option3"
             ),
             facetName = c(as.list(rep("Higher biomass scenario", 3)), as.list(rep("Lower biomass scenario", 3)), as.list(rep("Other biomass scenario", 3))),
             chooseArea = 0,
             percentile = c(0.025, 0.975),
             doHist = TRUE,
             dpi = 300)
