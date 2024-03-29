% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YPR.R
\name{gtgYPRWrapper}
\alias{gtgYPRWrapper}
\title{YPR eumetric simulation using growth-type group model}
\usage{
gtgYPRWrapper(
  LifeHistoryObj,
  LcStep = 1,
  F_MStep = 0.2,
  waitName = NULL,
  hostName = NULL,
  gtg = 13,
  stepsPerYear = 12,
  selType = NULL,
  selParams = NULL,
  inchStep = NULL
)
}
\arguments{
\item{LifeHistoryObj}{A life history object.}

\item{LcStep}{Length step size in cm for sequence of length at vulnerability. Approx. knife edge vul. SL50 = Lc, SL95 = Lc + 1}

\item{F_MStep}{F/M ratio step size for sequence of F_M}

\item{waitName}{When used within a shiny app, this function can update a host from the waiter package. See example.}

\item{hostName}{When used within a shiny app, this function can update a host from the waiter package. See example.}

\item{gtg}{The number of growth-type groups. Default is 13.}

\item{stepsPerYear}{The number of steps per year. Default is 12.}

\item{selType}{When selectivity specified type stated here.}

\item{selParams}{Parameters for selectivity function}

\item{inchStep}{Optional parameter adding Lc corresponding to inch-based increments}
}
\description{
Produces YPR, Yield, and SPR arrays across combinations of F_M and Lc
}
\details{
Utilizes the growth-type group model to calculate YPR, Yield, and SPR across factorial combination of F_M and Lc.
This produces an array of values for each output, which is useful for surface plot, finding MSY, etc.
Yield is presented in relative terms, with a maximum of 1. YPR is not presented in relative terms, but the user
may want to re-calculate this quantity as relative YPR for presentation purposes.
Required parameters of the LifeHistory object are: Linf, L50, L95delta, M, K, Steep, LW_A, LW_B
}
\examples{
lh<-new("LifeHistory")
lh@Linf<-100
lh@L50<-66
lh@L95delta<-1
lh@MK<-1.5
lh@LW_A<-0.01
lh@LW_B<-3
#' sim<-gtgYPRWrapper(lh)

#################################################
#Use of hostess loading bar from waiter package
################################################
\dontrun{
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


shinyApp(ui, server)}
}
