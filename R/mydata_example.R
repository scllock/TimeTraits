#' Arabidopsis thaliana luciferase circadian rhythm dataset
#'
#' Raw luminescence time-series data from an Arabidopsis thaliana luciferase
#' reporter experiment examining circadian responses to changes in photoperiod.
#'
#' Plants (Ws-2 wild-type and phyB mutant lines) were entrained under short-day
#' conditions (8 hours light: 16 hours dark) and subsequently shifted to long-day
#' conditions (16 hours light: 8 hours dark). Luminescence was recorded at regular
#' time intervals for approximately 166 hours.
#'
#' The dataset is provided in wide format, with one column per individual plant
#' and one row per time point (h).
#'
#' @docType data
#'
#' @usage data(mydata_example)
#'
#' @format A data frame with 300 rows and 189 variables:
#' \describe{
#'   \item{new_time}{Numeric vector giving time since experiment start (hours).}
#'   \item{WS2.001--WS2.096}{Luminescence counts per second for individual WS-2
#'     wild-type plants. Each column corresponds to one biological replicate.}
#'   \item{phyb.001--phyb.096}{Luminescence counts per second for individual phyB
#'     mutant plants. Each column corresponds to one biological replicate.}
#' }
#'
#' All luminescence values are raw photon counts per second (cps).
#'
#' @source Experimental data generated as part of a luciferase circadian rhythm study.
#'
#'
#' @keywords circadian luminescence time-series functional-data-analysis
#'
#' @examples
#' data(mydata_example)
#' head(mydata_example)
"mydata_example"

