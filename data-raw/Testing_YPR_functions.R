

devtools::load_all()

#----------------------
#Life history demo
#---------------------
lh<-LifeHistoryExample
lh@MK<-2
lh@M<-0.2
lh@K<-0.1
lh@Tmax<-1 #Keeping this set at 1 will override Tmax and use -log(0.01)/M for max age

ta<-new("TimeArea")
ta@gtg<-13
ta@stepsPerYear<-12

ptm<-proc.time()
x<-LHwrapper(LifeHistoryObj = lh, TimeAreaObj=ta)
print("Time in minutes: ")
print((proc.time()-ptm)/60)


#----------------------------------
#Demonstrate selectivity plotting
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

ta<-new("TimeArea")
ta@gtg<-13
ta@stepsPerYear<-1

X<-selWrapper(LifeHistoryObj=lh, TimeAreaObj = ta, FisheryObj = FisheryObj, doProjection = TRUE, doPlot = TRUE)

#----------------------------------------
#Demonstrate plotting function of solveD
#----------------------------------------
FisheryObj<-new("Fishery")
FisheryObj@title<-"Test"
FisheryObj@historicalVulType<-"logistic"
FisheryObj@historicalVulParams<-c(50,75)
FisheryObj@historicalRetType<-"full"
FisheryObj@historicalRetParams<-c(60,70)
FisheryObj@historicalRetMax<-1
FisheryObj@historicalDmort<-0

ta<-new("TimeArea")
ta@gtg<-13
ta@stepsPerYear<-12

lh<-LHwrapper(LifeHistoryObj=LifeHistoryExample, TimeAreaObj = ta)
sel<-selWrapper(LifeHistoryObj=LifeHistoryExample, TimeAreaObj = ta, FisheryObj = FisheryObj, doProjection = TRUE, doPlot = TRUE)

X<-solveD(lh, sel, doFit = FALSE, F_in = 0, doPlot = TRUE)




#---------------------
#Comparing YPR curves
#---------------------
#LH
lh<-LifeHistoryExample
lh@MK<-0.5
lh@M<-0.1
lh@K<-0.2
lh@Tmax<-1 #Keeping this set at 1 will override Tmax and use -log(0.01)/M for max age

ptm<-proc.time()
sim<-lbsprSimWrapper(LifeHistoryObj = lh)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

ptm<-proc.time()
sim2<-lbsprSimWrapperAbsel(LifeHistoryObj = lh)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

ptm<-proc.time()
sim3<-gtgYPRWrapper(LifeHistoryObj = lh, gtg=13, stepsPerYear = 1)
print("Time in minutes: ")
print((proc.time()-ptm)/60)

ptm<-proc.time()
sim4<-gtgYPRWrapper(LifeHistoryObj = lh, gtg=21, stepsPerYear = 12)
print("Time in minutes: ")
print((proc.time()-ptm)/60)


#------------
#F_M on x axis
#-------------
#YPR
plot(sim@sim$F_M, sim@sim$YPR_EU[,70]/max(sim@sim$YPR_EU), type="l", col="red")
lines(sim2@sim$F_M, sim2@sim$YPR_EU[,70]/max(sim2@sim$YPR_EU), col="blue")
lines(sim3@sim$F_M, sim3@sim$YPR_EU[,70]/max(sim3@sim$YPR_EU), col="green")
lines(sim4@sim$F_M, sim4@sim$YPR_EU[,70]/max(sim4@sim$YPR_EU), col="orange")

#Yield
plot(sim@sim$F_M, sim@sim$Yield_EU[,70], type="l", col="red", main = "0.01")
lines(sim2@sim$F_M, sim2@sim$Yield_EU[,70], col="blue")
lines(sim3@sim$F_M, sim3@sim$Yield_EU[,70]/max(sim3@sim$Yield_EU), col="green")
lines(sim4@sim$F_M, sim4@sim$Yield_EU[,70]/max(sim4@sim$Yield_EU), col="orange")

#SPR
plot(sim@sim$F_M, sim@sim$SPR_EU[,70], type="l", col="red", main = "0.01")
lines(sim2@sim$F_M, sim2@sim$SPR_EU[,70], col="blue")
lines(sim3@sim$F_M, sim3@sim$SPR_EU[,70], col="green")
lines(sim4@sim$F_M, sim4@sim$SPR_EU[,70], col="orange")

#------------
#Lc on x axis
#-------------
#YPR
plot(sim@sim$Lc, sim@sim$YPR_EU[15,]/max(sim@sim$YPR_EU), type="l", col="red", ylim = c(0,1))
lines(sim2@sim$Lc, sim2@sim$YPR_EU[15,]/max(sim2@sim$YPR_EU), col="blue")
lines(sim3@sim$Lc, sim3@sim$YPR_EU[15,]/max(sim3@sim$YPR_EU), col="green")
lines(sim4@sim$Lc, sim4@sim$YPR_EU[15,]/max(sim4@sim$YPR_EU), col="orange")

#Yield
plot(sim@sim$Lc, sim@sim$Yield_EU[15,], type="l", col="red", main = "0.01")
lines(sim2@sim$Lc, sim2@sim$Yield_EU[15,], col="blue")
lines(sim3@sim$Lc, sim3@sim$Yield_EU[15,]/max(sim3@sim$Yield_EU), col="green")
lines(sim4@sim$Lc, sim4@sim$Yield_EU[15,]/max(sim4@sim$Yield_EU), col="orange")

#SPR
plot(sim@sim$Lc, sim@sim$SPR_EU[15,], type="l", col="red", main = "0.01")
lines(sim2@sim$Lc, sim2@sim$SPR_EU[15,], col="blue")
lines(sim3@sim$Lc, sim3@sim$SPR_EU[15,], col="green")
lines(sim4@sim$Lc, sim4@sim$SPR_EU[15,], col="orange")


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
  gtgYPRWrapper(LifeHistoryObj = LifeHistoryExample, waitName=w, hostName=host)
  w$hide()

}


shinyApp(ui, server)
