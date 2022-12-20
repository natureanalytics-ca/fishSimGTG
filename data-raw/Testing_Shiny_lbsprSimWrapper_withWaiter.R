
library(shiny)
library(waiter)
library(fishSimGTG)
library(bs4Dash)
library(here)
devtools::load_all()

ui <- fluidPage(
  br(),
  actionButton(
    "rn",
    "Run",
    status = "danger"
  ),
  useWaiter(),
  useHostess() # include dependencies

)

server <- function(input, output){

  #--------------------------------
  #Loaders
  #--------------------------------
  host <- Hostess$new()
  waitLoad <- Waiter$new(
    html = host$get_loader(
      preset = "bubble",
      text_color = "black",
      center_page = TRUE,
      class = "",
      min = 0,
      max = 100,
      svg = NULL,
      progress_type = "fill",
      fill_direction = c("btt", "ttb", "ltr", "rtl"),
      stroke_direction = c("normal", "reverse"),
      fill_color = NULL,
      stroke_color = NULL
    ),
    color = transparent(alpha = 0.8),
    fadeout = TRUE
  )


  observeEvent(input$rn, {

    #-------------------
    #LBSPR example
    #-------------------

    #X<-runIt
    # sim<-fishSimGTG::lbsprSimWrapper(LifeHistory = LifeHistoryExample,
    #                                  binWidth = 2,
    #                                  binMin = 0,
    #                                  LcStep = 2,
    #                                  F_MStep = 0.2,
    #                                  waitName = waitLoad,
    #                                  hostName = host
    # )

    #-------------------
    #MSE sims example
    #-------------------


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
    TimeAreaObj@iterations = 50
    TimeAreaObj@historicalYears = 10
    TimeAreaObj@historicalBio = 0.5
    TimeAreaObj@historicalBioType = "relB"
    TimeAreaObj@move <- matrix(c(1,0, 0,1), nrow=2, ncol=2, byrow=FALSE)
    TimeAreaObj@historicalEffort<-matrix(1:1, nrow = 10, ncol = 2, byrow = FALSE)

    HistFisheryObj<-new("Fishery")
    HistFisheryObj@title<-"Example"
    HistFisheryObj@vulType<-"logistic"
    HistFisheryObj@vulParams<-c(40.1,0.1) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
    HistFisheryObj@retType<-"full"
    HistFisheryObj@retMax <- 1
    HistFisheryObj@Dmort <- 0

    StochasticObj<-new("Stochastic")
    StochasticObj@historicalBio = c(0.3, 0.6)

    ProFisheryObj<-new("Fishery")
    ProFisheryObj@title<-"Example"
    ProFisheryObj@vulType<-"logistic"
    ProFisheryObj@vulParams<-c(40.1,0.1) #Approx. knife edge based on input value of 40.1. Must put slightly higher value for second parameter
    ProFisheryObj@retType<-"full"
    #ProFisheryObj@retParams <- c(35.6, 0.1)
    ProFisheryObj@retMax <- 1
    ProFisheryObj@Dmort <- 0
    ProFisheryObj_list<-list(ProFisheryObj, ProFisheryObj)

    StrategyObj <- new("Strategy")
    StrategyObj@projectionYears <- 50
    StrategyObj@projectionName<-"projectionStrategy"
    StrategyObj@projectionParams<-list(bag = c(-99,-99), effort = matrix(1:1, nrow=50, ncol=2, byrow = FALSE), CPUE = c(5,8), CPUEtype = "retN", effortImpError = c(0.7, 1.3))

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
                  seed = 10,
                  waitName = waitLoad,
                  hostName = host
    )

  })

  runIt<-function(waitName = NULL, hostName=NULL){
    if(!is.null(hostName) & !is.null(waitName)){
      waitName$show()
    }
    for(i in 1:10){
      Sys.sleep(0.2) # random sleep
      if(!is.null(hostName) & !is.null(waitName)){
        hostName$set(i * 10)
      }
    }
    if(!is.null(hostName) & !is.null(waitName)){
      waitName$hide()
    }
    x<-10*2
    return(x)
  }

}

shinyApp(ui, server)
