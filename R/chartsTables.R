

#---------------------------------------
#Relative biomass vs. relative catch
#---------------------------------------

#Roxygen header
#'Relative change in biomass vs. relative catch scatter plot
#'
#'Plots the ratio of year prior to projection to a specified projection year. Iterations are summarized as the median ratio. Calculation are made for relative SSB (y-axis) and relative catch in weight (a-axis). Saves a chart to the wd. Examples can be found in data-raw/DAR_projection_example.R
#'
#' @param wd A working directly where the output of runProjection is saved
#' @param fileName List. List of file names in the working directory
#' @param facetName List. Plot uses ggplot2 facets. File names can be assigned to facets, thus producing separate plots for a given facet, such a low M scenario vs. a high M scenario.
#' @param chooseArea The area to summarized in the plot. Value of 0 sums quantities across all areas
#' @param proYear Numeric. The projection year used to calculate relative change in SSB and catch in weight.
#' @param dpi Resolution in dots per inch of the resulting saved chart.
#' @param imageName Character. A name for the resulting plot(s)
#' @import ggplot2 ggrepel
#' @export

relSSBscatter<-function(wd, fileName, facetName, chooseArea = 0, proYear, dpi = 300, imageName = "Relative_SSB_catch"){

  catchB_median<-SSB_median<-nm<-NULL

  totalSSB<-data.frame()
  areaSSB<-data.frame()
  for(i in 1:NROW(fileName)){
    data<-readRDS(paste0(wd, "/", fileName[[i]],".rds"))

    refYear<-1 + data$TimeAreaObj@historicalYears
    getYear<-1 + data$TimeAreaObj@historicalYears + proYear

    #Choose area
    if(chooseArea == 0) { #Sum across areas
      SSB_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        sum(data$dynamics$SB[getYear,x,]) / sum(data$dynamics$SB[refYear,x,])
      })
      cat_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        sum(data$dynamics$catchB[getYear,x,]) / sum(data$dynamics$catchB[refYear,x,])
      })
      totalSSB<-rbind(totalSSB, list(fct = facetName[[i]], nm = data$titleStrategy, SSB_median = median(SSB_tmp), catchB_median = median(cat_tmp)))
    } else {#Specific area selected by user
      SSB_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        data$dynamics$SB[getYear,x,chooseArea] / data$dynamics$SB[refYear,x,chooseArea]
      })
      cat_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        data$dynamics$catchB[getYear,x,chooseArea] / data$dynamics$catchB[refYear,x,chooseArea]
      })
      totalSSB<-rbind(totalSSB, list(fct = facetName[[i]], nm = data$titleStrategy, SSB_median = median(SSB_tmp), catchB_median = median(cat_tmp)))
    }
  }

  #Total relative SSB
  ggplot(totalSSB, aes(catchB_median, SSB_median)) +
    geom_vline(xintercept = 1, color="lightgrey") +
    geom_hline(yintercept = 1, colour="lightgrey") +
    geom_point(size=1) +
    geom_label_repel(data=totalSSB,
                    aes(catchB_median, SSB_median, label = nm),
                    fontface = 'bold',
                    color="cornflowerblue",
                    size=2,
                    segment.size = 0.1,
                    segment.color = "black",
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0, "lines"),
                    show.legend = FALSE,
                    max.overlaps = 25
    ) +
    ylab(paste("Relative biomass after", proYear, "years")) +
    xlab(paste("Relative catch weight after", proYear, "years")) +
    theme_classic() +
    theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
          strip.text.y=element_text(colour = "black", size=8, face="bold"),
          strip.background = element_rect(fill ="lightgrey"),
          axis.text=element_text(size=6),
          panel.border = element_rect(linetype = "solid", colour = "black", fill=NA))+
    facet_wrap(~fct, ncol=2)
  ggsave(filename = paste0(wd, "/", imageName, ".png"), device = "png", dpi = dpi, width = 6, height = (3*NROW(unique(totalSSB$fct))/2), units = "in")
}

#---------------------------------------
#Relative biomass time series
#---------------------------------------

#Roxygen header
#'Relative change in biomass vs. relative catch scatter plot
#'
#'Time series of changes in spawning stock biomass and catch weight (two separate plots produced) relative to the year prior to projection time period. Saves a chart to the wd. Examples can be found in data-raw/DAR_projection_example.R
#'
#' @param wd A working directly where the output of runProjection is saved
#' @param fileName List. List of file names in the working directory
#' @param facetName List. Plot uses ggplot2 facets. File names can be assigned to facets, thus producing separate plots for a given facet, such a low M scenario vs. a high M scenario.
#' @param chooseArea The area to summarized in the plot. Value of 0 sums quantities across all areas
#' @param percentile Vector of length two indicating centered percent of observations to represent uncertainty in time series plots. For example, 95% centered observations should be entered as c(0.025, 0.975)
#' @param doHist Should the historical time period be included in the plot?
#' @param dpi Resolution in dots per inch of the resulting saved chart.
#' @param imageName Character. A name for the resulting plot(s)
#' @import ggplot2 ggrepel
#' @importFrom stats quantile
#' @export

relSSBseries<-function(wd, fileName, facetName, chooseArea = 0, percentile = c(0.025, 0.975), doHist = FALSE, dpi = 300, imageName = "Relative_timeSeries"){

  year<-med<-lower<-upper<-NULL

  totalSSB<-data.frame()
  totalcatchB<-data.frame()
  for(i in 1:NROW(fileName)){

    data<-readRDS(paste0(wd, "/", fileName[[i]],".rds"))
    refYear <- 1 + data$TimeAreaObj@historicalYears
    startYear <- ifelse(doHist, 1, refYear)
    endYear <- 1 + data$TimeAreaObj@historicalYears + data$StrategyObj@projectionYears

    #Choose area
    if(chooseArea == 0) { #Sum across areas
      SSB_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        rowSums(data$dynamics$SB[startYear:endYear,x,])/sum(data$dynamics$SB[refYear,x,])
      })

      cat_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        rowSums(data$dynamics$catchB[startYear:endYear,x,])/sum(data$dynamics$catchB[refYear,x,])
      })
    } else { #Specific area selected by user
      SSB_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        data$dynamics$SB[startYear:endYear,x,chooseArea]/data$dynamics$SB[refYear,x,chooseArea]
      })

      cat_tmp<-sapply(1:data$TimeAreaObj@iterations, FUN=function(x){
        data$dynamics$catchB[startYear:endYear,x,chooseArea]/data$dynamics$catchB[refYear,x,chooseArea]
      })
    }

    SSB_range<-t(sapply(1:NROW(SSB_tmp), FUN=function(x){
      c(quantile(SSB_tmp[x,], probs =  percentile), median(SSB_tmp[x,]), x)
    }))
    totalSSB<-rbind(totalSSB, list(fct = rep(facetName[[i]], NROW(SSB_range)), nm = rep(data$titleStrategy, NROW(SSB_range)), lower =  SSB_range[,1], upper = SSB_range[,2], med = SSB_range[,3], year = SSB_range[,4]))

    cat_range<-t(sapply(1:NROW(cat_tmp), FUN=function(x){
      c(quantile(cat_tmp[x,], probs =  percentile), median(cat_tmp[x,]), x)
    }))
    totalcatchB<-rbind(totalcatchB, list(fct = rep(facetName[[i]], NROW(cat_range)), nm = rep(data$titleStrategy, NROW(cat_range)), lower =  cat_range[,1], upper = cat_range[,2], med = cat_range[,3], year = cat_range[,4]))
  }

  #Total relative SSB
  ggplot(totalSSB, aes(x = year, y = med)) +
    geom_line(color = "cornflowerblue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    ylab("Relative biomass") +
    xlab("Year") +
    theme_classic() +
    theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
          strip.text.y=element_text(colour = "black", size=8, face="bold"),
          strip.background = element_rect(fill ="lightgrey"),
          axis.text=element_text(size=6),
          panel.border = element_rect(linetype = "solid", colour = "black", fill=NA)) +
    facet_wrap(~ fct + nm, ncol=NROW(unique(totalSSB$nm)), scale = "free_y")
  ggsave(filename = paste0(wd, "/", imageName, "_SSB.png"), device = "png", dpi = dpi, width = 6, height = 3.5*NROW(unique(totalSSB$fct)), units = "in")

  #Total relative catch weight
  ggplot(totalcatchB, aes(x = year, y = med)) +
    geom_line(color = "cornflowerblue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    ylab("Relative catch weight") +
    xlab("Year") +
    theme_classic() +
    theme(strip.text.x=element_text(colour = "black", size=8, face="bold"),
          strip.text.y=element_text(colour = "black", size=8, face="bold"),
          strip.background = element_rect(fill ="lightgrey"),
          axis.text=element_text(size=6),
          panel.border = element_rect(linetype = "solid", colour = "black", fill=NA)) +
    facet_wrap(~ fct + nm, ncol=NROW(unique(totalcatchB$nm)), scale = "free_y")
  ggsave(filename = paste0(wd, "/", imageName, "_catchB.png"), device = "png", dpi = dpi, width = 6, height = 3.5*NROW(unique(totalcatchB$fct)), units = "in")
}
