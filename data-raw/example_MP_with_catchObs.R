rm(list=ls())
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

    #Need to combine obs (e.g., index and catch) in a list before returning
    combined_data <- list()

    # get index obs data and add to combined_data list
    if(!is.null(IndexObj)) {
      index_result <- calculate_single_Index(dataObject)
      # add all index columns
      for(col_name in names(index_result)) {
        combined_data[[col_name]] <- index_result[[col_name]]
      }
    }

    # get catch obs data and add to combined_data list
    if(!is.null(CatchObsObj)) {
      catch_result <- calculate_single_CatchObs(dataObject)
      # add all catch columns
      for(col_name in names(catch_result)) {
        combined_data[[col_name]] <- catch_result[[col_name]]
      }
    }

    #User defines variable names for returned list, as this info will be used in Phase 2 (HCR)
    #return(list()) #e.g., return(list(year = j, iteration = k, CPUE = 1.4))
    return(combined_data)
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

#-----------------------------------------------#
#------------NEW--------------------------------#
#-----------------------------------------------#
# Create a new Catchobs class
catch_obs1 <- new("CatchObs",
                      catchID = "Fishery Catch",
                      title = "Fishery catch observations",
                      areas = c(1, 2),  # Both areas like CPUE
                      catchYears = 2:21,  # Same years as CPUE
                      reporting_rates = c(seq(0.6, 1, length=10), rep(1, 10)),  # Improving historical, non systematic bias in projection
                      obs_CVs = matrix(c(rep(0.2, 20), rep(0.4, 20)), ncol=2))  # CV bounds: 0.2-0.4 for all years


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB1,
  CatchObsObj = catch_obs1,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_with_catch_obs",
  doPlot = TRUE
)

out1<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_with_catch_obs")
out1$HCR$decisionData
out1$HCR$decisionData$CPUE_1
out1$HCR$decisionData$observed_catch
out1$dynamics$Ftotal


#--------------------------------------#
# Plot function for observation models #
#--------------------------------------#

#testing plot fucntion
#tibble_data=X1$HCR$decisionData


plotIndex_tibble <- function(tibble_data, save_plot = FALSE,
                             filename = "index_plot.jpeg") {

  # extract data from the tibble
  title <- unique(tibble_data$title)   #unique gets unique values (removes duplicates)
  historical_end <- unique(tibble_data$historical_end)

  # identify CPUE/Survey columns (the actual index values)
  # find the columns names
  # find the pattern name with grep (must start (^) with CPUE_|Survey_    and have one or more digists(\\d+) $= ends here)
  index_columns <- grep("^(CPUE_|Survey_)\\d+$", names(tibble_data), value = TRUE)

  # create a empty data frame to prepare data for plotting
  all_data <- data.frame()


  # index_col (a column) = eg. "CPUE_1"
  for(index_col in index_columns) {
    # get the corresponding columns
    indextype_col <- paste0(index_col, "_indextype")
    areas_col <- paste0(index_col, "_areas")
    indexyears_col <- paste0(index_col, "_indexYears")

    # extract indexYears for this index
    index_years_str <- unique(tibble_data[[indexyears_col]])[1]
    #  "2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21"
    index_years <- as.numeric(unlist(strsplit(index_years_str, "_")))  # Splits string by "_"
    #2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    # get index type and areas for panel naming
    index_type <- unique(tibble_data[[indextype_col]])[1]
    areas_str <- unique(tibble_data[[areas_col]])[1]

    # create panel name
    panel_name <- paste(index_col, "-", "Area(s)", areas_str)

    # prepare data for this index
    index_data <- tibble_data %>%
      select(k, j, all_of(index_col)) %>%
      rename(
        iteration = k,
        year = j,
        value = !!sym(index_col) #!!sym=creates a "symbol" representing the column name CPUE_1- !! =unquote = CPUE_1
      ) %>%
      #filter(!is.na(value)) %>%
      mutate(
        panel = panel_name,
        type = "Iterations",
        iteration_label = paste0("Iter_", iteration)
      )

    # calculate median across iterations for each year
    # drops remove the grouping after summarising
    median_data <- index_data %>%
      group_by(year, panel) %>%
      summarise(
        value = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      #text labels created using later for plotting
      mutate(
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )

    # set NA for years not in indexYears (to match original behavior)
    all_years <- min(tibble_data$j):max(tibble_data$j)
    missing_years <- setdiff(all_years, index_years)

    if(length(missing_years) > 0) {
      missing_data <- data.frame(
        year = missing_years,
        value = NA,
        panel = panel_name,
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )
      median_data <- rbind(median_data, missing_data)
    }

    # combine iteration and median data
    combined_data <- rbind(
      index_data %>% select(year, value, panel, type, iteration_label),
      median_data %>% select(year, value, panel, type, iteration_label)
    )

    all_data <- rbind(all_data, combined_data)
  }
  #Y-axis label based on data type
  use_weight <- unique(tibble_data$useWeight)[1]
  final_indextype <- unique(tibble_data$final_indextype)[1]

  if(final_indextype == "FD") {
    # Fishery Dependent (CPUE)
    y_label <- if(use_weight) "CPUE (biomass)" else "CPUE (numbers)"
  } else if(final_indextype == "FI") {
    # Fishery Independent (Survey)
    y_label <- if(use_weight) "Survey (biomass)" else "Survey (numbers)"
  } else {
    # Mixed types
    y_label <- if(use_weight) "Index (biomass)" else "Index (numbers)"
  }



  # create the plot (matching original style)
  p <- ggplot(all_data, aes(x = year, y = value)) +
    geom_line(data = subset(all_data, type == "Iterations"),
              aes(group = iteration_label),
              color = "steelblue", alpha = 0.6, size = 0.5) +
    geom_point(data = subset(all_data, type == "Iterations" & !is.na(value)),
               color = "steelblue", alpha = 0.7, size = 1) +
    geom_line(data = subset(all_data, type == "Median"),
              color = "black", size = 1.2) +
    geom_point(data = subset(all_data, type == "Median" & !is.na(value)),
               color = "black", size = 1.5) +
    facet_wrap(~ panel, scales = "free_y") +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    labs(title = title, x = "Year", y = y_label) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )


  if(save_plot) {
    ggsave(filename, plot = p, width = 12, height = 8, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}



#--------------------------------------------#
# Plot function for catch observation models #
#--------------------------------------------#

plotCatch_tibble <- function(tibble_data, save_plot = FALSE,
                             filename = "catch_plot.jpeg") {

  # extract data from the tibble
  title <- unique(tibble_data$catchID)[1]
  historical_end <- unique(tibble_data$end_hist)

  # only plot observed_catch
  index_columns <- "observed_catch"

  # create a empty data frame
  all_data <- data.frame()

  # process observed_catch only
  for(index_col in index_columns) {

    # define panel name
    areas_str <- unique(tibble_data$areas_included)[1]
    panel_name <- paste("Observed catch -", "Area(s)", areas_str)

    # years with catch data
    index_years <- sort(unique(tibble_data$j[!is.na(tibble_data[[index_col]])]))

    # prepare data same structure as index plot
    index_data <- tibble_data %>%
      select(k, j, all_of(index_col)) %>%
      rename(
        iteration = k,
        year = j,
        value = !!sym(index_col)
      ) %>%
      mutate(
        panel = panel_name,
        type = "Iterations",
        iteration_label = paste0("Iter_", iteration)
      )

    # calculate median across iterations for each year as index plot
    median_data <- index_data %>%
      group_by(year, panel) %>%
      summarise(
        value = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )

    # set NA for years not in catchYears (as index plot)
    all_years <- min(tibble_data$j):max(tibble_data$j)
    missing_years <- setdiff(all_years, index_years)

    if(length(missing_years) > 0) {
      missing_data <- data.frame(
        year = missing_years,
        value = NA,
        panel = panel_name,
        type = "Median",
        iteration = "Median",
        iteration_label = "Median"
      )
      median_data <- rbind(median_data, missing_data)
    }

    # combine iteration and median data (as index plot)
    combined_data <- rbind(
      index_data %>% select(year, value, panel, type, iteration_label),
      median_data %>% select(year, value, panel, type, iteration_label)
    )

    all_data <- rbind(all_data, combined_data)
  }

  # Y-axis label for catch
  y_label <- "Observed Catch (biomass)"

  # create the plot (as index plot)
  p <- ggplot(all_data, aes(x = year, y = value)) +
    geom_line(data = subset(all_data, type == "Iterations"),
              aes(group = iteration_label),
              color = "steelblue", alpha = 0.6, size = 0.5) +
    geom_point(data = subset(all_data, type == "Iterations" & !is.na(value)),
               color = "steelblue", alpha = 0.7, size = 1) +
    geom_line(data = subset(all_data, type == "Median"),
              color = "black", size = 1.2) +
    geom_point(data = subset(all_data, type == "Median" & !is.na(value)),
               color = "black", size = 1.5) +
    facet_wrap(~ panel, scales = "free_y") +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    labs(title = title, x = "Year", y = y_label) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if(save_plot) {
    ggsave(filename, plot = p, width = 12, height = 8, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}




p1 <- plotIndex_tibble(out1$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEB1.jpeg")

p1


p2 <- plotCatch_tibble (out1$HCR$decisionData,
                           save_plot = TRUE,
                             filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/plot_catchobs1.jpeg")
p2

#higer F in the projection period exaplain the increases in catchs
out1$dynamics$Ftotal


#---------------------------------------------------------------------------------------#
# EXAMPLE 2                                                                             #
# Observation model 2: Simulation of two CPUE (biomass) index covering both areas       #
# adding catch observtion model                                                         #
#---------------------------------------------------------------------------------------#

#Define the new Index class
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


# Define CVs for CatchObs class
total_years <- strategy@projectionYears+ta@historicalYears+1
end_hist <- ta@historicalYears+1
cv_periods <- matrix(nrow = total_years, ncol = 2)

for(j in 2:total_years) {
  if(j <= end_hist) {
    cv_periods[j, ] <- c(0.2, 0.4)  # Historical: CV between 0.2-0.4
  } else {
    cv_periods[j, ] <- c(0.1, 0.2)  # Projection: CV between 0.1-0.2
  }
}

#Define the new CatchObs class
catch_obs2 <- new("CatchObs",
                  catchID = "period_specific_CV",
                  title = "Different CV bounds by period",
                  areas = c(1, 2),
                  catchYears = 2:total_years,
                  #improving reporting rates over time
                  reporting_rates = seq(0.5, 0.9, length.out = total_years-1),
                  #more precise catch record for the projection period
                  obs_CVs = cv_periods)


# Now run the simulation
runProjection(
  LifeHistoryObj = lh,
  TimeAreaObj = ta,
  HistFisheryObj = fishery,
  ProFisheryObj_list = list(ProFisheryObj, ProFisheryObj),
  StrategyObj = strategy,
  StochasticObj = stochastic,
  IndexObj = cpueB2,
  CatchObsObj = catch_obs2,
  customToCluster = "testMP",
  wd = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole",
  fileName = "test_run_MP_with_catch_obs2",
  doPlot = TRUE
)

out2<-readProjection("P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole", "test_run_MP_with_catch_obs2")
out2$HCR$decisionData
out2$HCR$decisionData$CPUE_1
out2$HCR$decisionData$observed_catch
out2$dynamics$Ftotal


p2 <- plotIndex_tibble(out2$HCR$decisionData,
                       save_plot = TRUE,
                       filename = "P:/Nature_Analytics_work/Simulation_obs_models1/data-test/Kole/index_plot_CPUEB2.jpeg")

p2
