

devtools::load_all()

#---------------------
#Comparing YPR curves
#---------------------
#LH
lh<-LifeHistoryExample
lh@MK<-0.5
lh@M<-0.1
lh@K<-0.2
lh@Tmax<-1 #Keeping this set at 1 will override Tmax and use -log(0.01)/M for max age
lh@isHermaph<-FALSE

ptm<-proc.time()
sim<-lbsprSimWrapper(LifeHistoryObj = lh)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

ptm<-proc.time()
sim2<-lbsprSimWrapperAbsel(LifeHistoryObj = lh)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

ptm<-proc.time()
sim3<-gtgSimWrapper(LifeHistoryObj = lh)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

#YPR
plot(sim@sim$F_M, sim@sim$YPR_EU[,20]/max(sim@sim$YPR_EU), type="l", col="red")
lines(sim2@sim$F_M, sim2@sim$YPR_EU[,20]/max(sim2@sim$YPR_EU), col="blue")
lines(sim3@sim$F_M, sim3@sim$YPR_EU[,20]/max(sim3@sim$YPR_EU), col="green")

#Yield
plot(sim@sim$F_M, sim@sim$Yield_EU[,20], type="l", col="red", main = "0.01")
lines(sim2@sim$F_M, sim2@sim$Yield_EU[,20], col="blue")
lines(sim3@sim$F_M, sim3@sim$Yield_EU[,20]/max(sim3@sim$Yield_EU), col="green")

#SPR
plot(sim@sim$F_M, sim@sim$SPR_EU[,69], type="l", col="red", main = "0.01")
lines(sim2@sim$F_M, sim2@sim$SPR_EU[,69], col="blue")
lines(sim3@sim$F_M, sim3@sim$SPR_EU[,69], col="green")


#----------------------------------------
#Demonstrate plotting function of solveD
#----------------------------------------
lh<-LifeHistoryExample
lh@MK<-2
lh@M<-0.2
lh@K<-0.1
lh@Tmax<-1 #Keeping this set at 1 will override Tmax and use -log(0.01)/M for max age

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@historicalVulType<-"logistic"
FisheryObj@historicalVulParams<-c(50,75)
FisheryObj@historicalRetType<-"full"
FisheryObj@historicalRetParams<-c(60,70)
FisheryObj@historicalRetMax<-1
FisheryObj@historicalDmort<-0

solveD(LifeHistoryObj=lh, FisheryObj=FisheryObj, doFit = FALSE, F_in = 0, doPlot = TRUE)

#----------------------------------
#Demonstrate selecvitity plotting
#-----------------------------------
lh<-LifeHistoryExample
lh@MK<-2
lh@M<-0.2
lh@K<-0.1
lh@Tmax<-1 #Keeping this set at 1 will override Tmax and use -log(0.01)/M for max age

FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@historicalVulType<-"logistic"
FisheryObj@historicalVulParams<-c(50,75)
FisheryObj@historicalRetType<-"slotLimit"
FisheryObj@historicalRetParams<-c(60,70)
FisheryObj@historicalRetMax<-1
FisheryObj@historicalDmort<-0.1
FisheryObj@projectionVulType<-"logistic"
FisheryObj@projectionVulParams<-c(50,75)
FisheryObj@projectionRetType<-"full"
FisheryObj@projectionRetMax<-0.5
FisheryObj@projectionDmort<-0

X<-selWrapper(LifeHistoryObj=lh, FisheryObj = FisheryObj, doProjection = TRUE, doPlot = TRUE)


#-----------------------
#Test spinner for Shiny
#-----------------------


library(shiny)
library(waiter)
library(fishSimGTG)

ui <- fluidPage(
 useWaiter(),
  useHostess(), # include dependencies
  actionButton("btn", "render"),
)

server <- function(input, output){

  host <- Hostess$new()
  w <- Waiter$new(
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
      stroke_color = "pink"
    ),
    color = transparent(alpha = 0.2),
    fadeout = TRUE
  )

 w$show()
  gtgSimWrapper(LifeHistoryObj = LifeHistoryExample, waitName=w, hostName=host)
  w$hide()

}


shinyApp(ui, server)
