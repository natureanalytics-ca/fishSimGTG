

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
#' @import tidyverse
#' @importFrom methods slot slotNames
#' @importFrom stats setNames
#' @export
#' @examples
#' library(tidyverse)
#' LH_to_dataframe(LifeHistoryExample)

LH_to_dataframe <- function(LifeHistoryObj, digits=3) {
  nms <- slotNames(LifeHistoryObj)
  lst <- lapply(nms, function(nm) slot(LifeHistoryObj, nm))
  ind<-which(lengths(lst)!=0)
  data.frame(setNames(lst[ind], nms[ind])) %>%
    mutate(across(c(Linf, K, t0, L50, L95, M, MK, LW_B, Steep), round, digits))  %>%
    mutate(across(c(Tmax), round, 1))  %>%
    mutate(across(c(LW_A), signif, digits))  %>%
    gather() %>%
    mutate_at("key", str_replace, "title", "Title") %>%
    mutate_at("key", str_replace, "speciesName", "Species") %>%
    mutate_at("key", str_replace, "shortDescription", "Short description") %>%
    mutate_at("key", str_replace, "L_type", "Length type") %>%
    mutate_at("key", str_replace, "L_units", "Length units") %>%
    mutate_at("key", str_replace, "Walpha_units", "Weight units") %>%
    mutate_at("key", str_replace, "Linf", "von Bertalanffy Loo") %>%
    mutate_at("key", str_replace, "K", "von Bertalanffy K") %>%
    mutate_at("key", str_replace, "t0", "von Bertalanffy t0") %>%
    mutate_at("key", str_replace, "L50", "Length at 50% maturity") %>%
    mutate_at("key", str_replace, "L95", "Length at 95% maturity") %>%
    mutate_at("key", str_replace, "M", "Natural mortality") %>%
    mutate_at("key", str_replace, "MK", "M/K") %>%
    mutate_at("key", str_replace, "LW_A", "Length-weight alpha") %>%
    mutate_at("key", str_replace, "LW_B", "Length-weight beta") %>%
    mutate_at("key", str_replace, "Tmax", "Maximum age") %>%
    mutate_at("key", str_replace, "Steep", "Beverton-Holt steepness") %>%
    mutate_at("key", str_replace, "author", "Author") %>%
    mutate_at("key", str_replace, "authAffiliation", "Author affiliation") %>%
    mutate_at("key", str_replace, "longDescription", "Long description")
}



