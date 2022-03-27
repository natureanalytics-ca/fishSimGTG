

#---------------------------------
#Life History object to data frame
#----------------------------------

#Roxygen header
#'Life History object to data frame
#'
#'Converts S4 life history object to data frame with nice variable names
#'
#' @param LifeHistoryObj  A life history object.
#' @param digits Integer indicating the number of decimal places (round). Applied as significant digits (signif) to L-W alpha.
#' @importFrom methods slot slotNames slot<-
#' @importFrom stats setNames
#' @importFrom tidyr gather
#' @import tidyverse
#' @export
#' @examples
#' library(tidyverse)
#' LH_to_dataframe(LifeHistoryExample)

LH_to_dataframe <- function(LifeHistoryObj, digits=3) {
  slot(LifeHistoryObj, 'appBuild', check=FALSE) <- NULL
  nms <- slotNames(LifeHistoryObj)
  lst <- lapply(nms, function(nm) slot(LifeHistoryObj, nm))
  ind<-which(lengths(lst)!=0)
  data.frame(setNames(lst[ind], nms[ind])) %>%
    mutate(across(any_of(c("Linf", "K", "t0", "L50", "L95", "M", "MK", "LW_B", "Steep", "H50", "H95", "recSD", "recRho")), round, digits))  %>%
    mutate(across(any_of(c("Tmax", "R0")), round, 1))  %>%
    mutate(across(any_of(c("LW_A")), signif, digits))  %>%
    rename_with(
      ~ case_when(
        . == "title" ~ "Title",
        . == "speciesName" ~ "Species",
        . == "shortDescription" ~ "Short description",
        . == "L_type" ~ "Length type",
        . == "L_units" ~ "Length units",
        . == "Walpha_units" ~ "Weight units",
        . == "Linf" ~ "von Bertalanffy Loo",
        . == "K" ~ "von Bertalanffy K",
        . == "t0" ~ "von Bertalanffy t0",
        . == "L50" ~ "Length at 50% maturity",
        . == "L95" ~ "Length at 95% maturity",
        . == "M" ~ "Natural mortality",
        . == "MK" ~ "M/K",
        . == "LW_A" ~ "Length-weight alpha",
        . == "LW_B" ~ "Length-weight beta",
        . == "Tmax" ~ "Maximum age",
        . == "Steep" ~ "Beverton-Holt steepness",
        . == "R0" ~ "Unfished recruitment",
        . == "recSD" ~ "Recruitment log-scale standard deviation",
        . == "recRho" ~ "Recruitment inter-annual correlation",
        . == "isHermaph" ~ "Protogynous hermaphrodite",
        . == "H50" ~ "Length at which cohort is 50% male",
        . == "H95" ~ "Length at which cohort is 95% male",
        . == "author" ~ "Author",
        . == "authAffiliation" ~ "Author affiliation",
        . == "longDescription" ~ "Long description",
        TRUE ~ .
      )
    ) %>%
    gather()
}



