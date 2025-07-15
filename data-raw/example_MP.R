# This first sections create the objects required by fishSimGTG and
# run the simulation using runProjection()

#clean the memory
rm(list=ls())
#load the modified fucntion
#source("obs_models_fishSimGTG_updated.R")
#load fishSimGTG
#library(fishSimGTG)
# used to test the fucntions to see if everything is working
devtools::load_all()
library(ggplot2)

# Create simple examples of each class to understand their structure
lh <- new("LifeHistory")
ta <- new("TimeArea")
fishery <- new("Fishery")
strategy <- new("Strategy")
stochastic <- new("Stochastic")

lh@title<-"Kole"
lh@speciesName<-"Ctenochaetus strigosus"
lh@Linf<-17.7
lh@K<-0.423
lh@t0<- -0.51
lh@L50<-8.4
lh@L95delta<-1.26
lh@M<-0.08
lh@L_type<-"FL"
lh@L_units<-"cm"
lh@LW_A<-0.046
lh@LW_B<-2.85
lh@Steep<-0.54
lh@recSD<-0 #Run with no rec var'n to see deterministic trends
lh@recRho<-0
lh@isHermaph<-FALSE
lh@R0<-10000

#Fishery ste up
fishery@title<-"Test"
fishery@vulType<-"logistic"
fishery@vulParams<-c(10.2,2) #Approx. knife edge
fishery@retType<-"full"
fishery@retMax <- 1
fishery@Dmort <- 0

#Time area set up
# The TimeArea class controls the size and structure of the simulation
ta@title = "Test"
ta@gtg = 13
ta@areas = 2
ta@recArea = c(0.99, 0.01)
ta@iterations = 2
ta@historicalYears = 10
ta@historicalBio = 0.5
ta@historicalBioType = "relB" # or SPR
ta@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
#Matrix of fishing effort multipliers
#Dimensions: [historicalYears × areas]
#Values are multipliers of equilibrium fishing rate
#e.g., 1.5 = 50% more fishing than equilibrium
ta@historicalEffort <- matrix(c(1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6),
                              nrow = 10, ncol = 2, byrow = FALSE)
dim(ta@historicalEffort)

#stochastic object set up
stochastic@historicalBio = c(0.3, 0.6)
stochastic@Steep= c(0.45, 0.75)

#projection fishery object
ProFisheryObj<-new("Fishery")
slotNames(ProFisheryObj)
ProFisheryObj@title<-"Test"
ProFisheryObj@vulType<-"logistic"
ProFisheryObj@vulParams<-c(10.2,0.1)
ProFisheryObj@retType<-"logistic"
ProFisheryObj@retParams <- c(10.2, 0.1)
ProFisheryObj@retMax <- 1
ProFisheryObj@Dmort <- 0



testMP <- function(phase, dataObject) {

  #------------------------------------
  #Unpacking: required - do not delete
  #------------------------------------
  #Unpack dataObject - objects passed from the OM to the MP
  for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])


  #------------------------------------------
  #Book keeping - be aware of year j indexing!
  #------------------------------------------
  #OM years are indexed 1:(1 + TimeAreaObj@historicalYears + StrategyObj@projectionYears)

  #Indexes specified within Class Strategy are typically indexed 1:projectionYears
  #To obtain the current projection year:
  #yr <- j - TimeAreaObj@historicalYears - 1

  #To obtain the terminal year of the historical period
  #yrHist <- TimeAreaObj@historicalYears + 1


  #-----------------------------------------------------------------
  #Summary of all available objects unpacked and passed from the OM
  #-----------------------------------------------------------------
  #j, current time step in OM
  #k, current iteration in OM
  #is, initial stock (equilibrium) for iteration k, year 1, contains values returned by solveD
  #lh, LHwrapper object for iteration k, contains values returned by LHwrapper
  #areas, contains TimeAreaObj@areas
  #ageClasses, contains the number of age classes. Obtained from lh$ageClasses. Caution: age 0 is always the first age class.
  #N, abundance-at-age for iteration k. A list of length lh$gtg that contains elements array(dim=c(ageClasses, years, areas))
  #selGroup, selectivity for iteration k and step j. A list of length areas with elements contains values returned by selWrapper
  #selHist, selectivity for iteration k for historical time period. A list of length areas with elements contains values returned by selWrapper
  #selPro, selectivity for iteration k for projection time period. A list of length areas with elements contains values returned by selWrapper
  #SB, spawing biomass at the beginnning of the year, array(dim=c(years, iterations, areas))
  #VB, vulnerable biomass at the beginning of the year, array(dim=c(years, iterations, areas))
  #RB, retained biomass at the beginning of the year, array(dim=c(years, iterations, areas))
  #catchN, annual catch in numbers in aggregate across gtg, array(dim=c(years, iterations, areas))
  #catchB, annual catch in biomass in aggregate across gtg, array(dim=c(years, iterations, areas))
  #discN, dead discards in numbers in aggregate across gtg, array(dim=c(years, iterations, areas))
  #discB, dead discards in biomass in aggregate across gtg, array(dim=c(years, iterations, areas))
  #Ftotal, fishing mortality rate, array(dim=c(years, iterations, areas))
  #SPR, spawning potential ratio, for entire population, array(dim=c(years, iterations))
  #relB, relative spawning biomass (depletion) for entire population, array(dim=c(years, iterations))
  #recN, annual age-0 recruitment in numbers, array(dim=c(years, iterations))
  #decisionData, optional data frame containing annual sampling details. (see phase 1 below)
  #decisionAnnual, optional data frame containing analysis and HCR info (see phase 2 below)
  #decisionLocal, mandatory data frame containing fishing mortality by area for year j and iteration k, (see phase 3 below)
  #RdevMatrix, recruitment deviations, contains values returned by recDev
  #Ddev, initial relative biomass deviations, contains values returned by bioDev
  #TimeAreaObj
  #StrategyObj


  #--------------------------------------
  #Phases of an MP
  #--------------------------------------
  #An MP can have up to 3 phases. Phase 3 is mandatory, phases 1 & 2 are optional
  #The OM will call this MP 3 times in each time step:
  #Phase 1: Observation model or simulated sampling/observations of the OM.
  #Phase 2: Analysis and harvest control rule
  #Phase 3: Calculate fishing mortality to be passed back to the OM.

  ########
  #Phase 1 - OPTIONAL. Observation process in year j
  #Phase 1 is called by the OM once for each iteration, k, and year, j.
  #Used to define imperfect sampling in year j
  #Returned list is appended to data frame: decisionData
  if(phase==1){
    #User defines computations needed.



    #User defines variable names for returned list, as this info will be used in Phase 2 (HCR)
    #return(list()) #e.g., return(list(year = j, iteration = k, CPUE = 1.4))
    return(calculate_single_Index(dataObject))
  }

  ########
  #Phase 2 - OPTIONAL. Analysis, decision-rule and/or HCR
  #Phase 2 is called by the OM once for each iteration, k, and year, j.
  #Used to define how catches will be regulated
  #Returned list is appended to data frame: decisionAnnual
  if(phase==2){
    #User defines computations needed (likely will rely on info saved to decisionData)

    #User defines variable names for returned list, as this info will be used in Phase 3
    return(list()) #e.g., return(list(year = j, iteration = k, TAC = 1400))

  }

  ########
  #Phase 3 - MANDATORY. Calculation of fishing mortality to be used by the OM
  #Phase 3 is called by the OM once for each iteration, k, and year, j.
  #Used to compute fishing mortality for the called k and j.
  #Returned list is appended to data frame: decisionLocal and retrieved by OM to calcualte catchN, catchB, etc.
  if(phase==3){

    #Must contain these elements, each a vector of length areas:
    year = rep(j, areas) #Time step to which fishing mortality applies
    iteration = rep(k, areas) #Iteration to which fishing mortality applies
    area = 1:areas #Area to which fishing mortality applies
    Flocal = rep(0.2, areas) #Calculation of fishing mortality (typically a complex computation is needed, e.g., to covert TAC to F)

    #Must return this list structure:
    return(list(year=year, iteration=iteration, area=area,  Flocal=Flocal))
  }
}


strategy@title <- "testing Obs"           # add a descriptive name
strategy@projectionYears <- 10                  #future projection
strategy@projectionName <- "testMP" # Which function to use - next: the parameters of projectionStrategy

#---------------------------------------------------------------------------------------#
# Observation model 1: Simulation of one CPUE (biomass) index covering both areas       #
#---------------------------------------------------------------------------------------#

# Create a new Index class
cpueB1 <- new("Index")
cpueB1@indexID <- "One cpueB1 - all_areas"
cpueB1@title <- "One cpueB1 - all_areas"
cpueB1@useWeight <- TRUE                # CPUE in biomass

# Define the survey design with cpue specific parameters
cpueB1@survey_design <- list(
  list(
    indextype = "FD",                          # new: add indextype to EACH index design
    areas = c(1, 2),                           # this CPUE covers both areas (1, 2)
    indexYears = 2:21,                         # all years have data (collect data every year)
    q_hist_bounds = c(0.00008, 0.00012),       # historical period - qbounds: 0.00008-0.00012
    q_proj_bounds = c(0.0001, 0.0003),         # projection period - qbounds: 0.0001-0.0003
    hyperstability_hist_bounds = c(0.85, 1.05),# historical - hypersbounds: 0.85-1.05
    hyperstability_proj_bounds = c(0.85, 1.05),# projection - hypersbounds: 0.85-1.05
    obsError_CV_hist_bounds = c(0.25, 0.35),   # historical CV bounds: 0.25-0.35
    obsError_CV_proj_bounds = c(0.2, 0.3)      # projection CV bounds: 0.2-0.3 (improved precision)
  )
)

# FD indices (CPUE) - empty list for FD data
cpueB1@selectivity_hist_list <- list()
cpueB1@selectivity_proj_list <- list()


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB1,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue1",
  doPlot = TRUE
)

X1<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue1")



#---------------------------------------------------------------------------------------#
# Observation model 2: Simulation of one CPUE (biomass) index covering both areas       #
#---------------------------------------------------------------------------------------#
cpueB2 <- new("Index")
cpueB2@indexID <- "Two_CPUE_programs_all_areas"
cpueB2@title <- "Two CPUE Programs - Both Cover All Areas"
cpueB2@useWeight <- TRUE                # CPUE in biomass

cpueB2@survey_design <- list(

  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1, 2),                              # covers all areas
    indexYears = 2:21,                            # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),           # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),           # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),     # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),   # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),        # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)       # projection CV: 0.15-0.25 (better precision)
  ),

  # CPUE 2
  list(
    indextype = "FD",
    areas = c(1, 2),                                     # also covers all areas
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    q_hist_bounds = c(0.0001, 0.00015),                  # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),                  # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),            # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),          # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),               # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)              # projection CV: 0.15-0.25 (some improvement)
  )
)

# FD indices (CPUE) - empty selectivity lists for FD data
cpueB2@selectivity_hist_list <- list()
cpueB2@selectivity_proj_list <- list()

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB2,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue2",
  doPlot = TRUE
)

X2<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue2")


#-----------------------------------------------------------------------
# Observation model 3: Two CPUE Biomass - Only Area 1
#-----------------------------------------------------------------------
cpueB3 <- new("Index")
cpueB3@indexID <- "Two_CPUE_programs_only_area1"
cpueB3@title <- "Two CPUE Programs - Cover only Area 1"
cpueB3@useWeight <- TRUE                # CPUE in biomass
cpueB3@survey_design <- list(
  # CPUE 1
  list(
    indextype = "FD",
    areas = c(1),                              # covers  area 1
    indexYears = 2:21,                      # annual data collection
    q_hist_bounds = c(0.0001, 0.00015),        # historical: 0.0001-0.00015
    q_proj_bounds = c(0.00015, 0.0003),        # projection: 0.00015-0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),  # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),# projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),     # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)    # projection CV: 0.15-0.25 (better precision)
  ),

  # CPUE 2
  list(
    indextype = "FD",
    areas = c(1),                                        # also covers Area 1
    indexYears = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20),  # data collected every other year
    q_hist_bounds = c(0.0001, 0.00015),                  # historical: 0.0001- 0.00015
    q_proj_bounds = c(0.00015, 0.0003),                  # projection: 0.00015- 0.0003
    hyperstability_hist_bounds = c(0.9, 1.1),            # historical: 0.9-1.1
    hyperstability_proj_bounds = c(0.85, 1.05),          # projection: 0.85-1.05
    obsError_CV_hist_bounds = c(0.2, 0.3),               # historical CV: 0.2-0.3
    obsError_CV_proj_bounds = c(0.15, 0.25)              # projection CV: 0.15-0.25 (some improvement)
  )
)
cpueB3@selectivity_hist_list <- list()
cpueB3@selectivity_proj_list <- list()

# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB3,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_Bcpue3",
  doPlot = TRUE
)

X3<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_Bcpue3")


