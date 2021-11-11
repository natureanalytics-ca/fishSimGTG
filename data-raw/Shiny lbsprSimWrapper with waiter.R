
library(shiny)
library(waiter)
library(fishSimGTG)
library(bs4Dash)

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

    #X<-runIt
    sim<-lbsprSimWrapper(LifeHistory = LifeHistoryExample,
                                     binWidth = 2,
                                     binMin = 0,
                                     LcStep = 2,
                                     F_MStep = 0.2,
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
