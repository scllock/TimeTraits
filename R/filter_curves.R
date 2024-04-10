#' Filter curves based on criteria
#'
#' The filter_curves function begins determining the time range for analysis and subsetting the data accordingly (default is full length of input). Samples with luminescence values below a specified threshold are filtered out (default is 0). The function applies linear regression to estimate and then remove trends from the curves, focusing on oscillatory components. Lomb-Scargle periodogram analysis is performed to detect rhythmic patterns in the detrended data (default is full length of input), followed by peak detection to identify significant periodic components. Samples failing to meet a significance threshold p > 0.001 in the periodogram analysis are filtered out. Finally, the function constructs and returns a filtered data matrix containing the specified time vector and curves that pass the filtering criteria, ready for further analysis or visualisation. Overall, \code{filter_curves} provides a comprehensive tool for preprocessing and analysing luminescence data, ensuring that only relevant and significant samples are retained for downstream analyses.
#'
#' @param x Matrix containing data values, each row is a time point and each column corresponds to a sample, can be a sparse matrix.
#' @param time Vector containing the time values.
#' @param from An integer specifying the time point to begin analysis. The default is the minimum of the time vector.
#' @param to An integer specifying the time point to end analysis. The default is the maximum of the time vector.
#' @param min_lum An integer specifying the minimum luminescence value to exclude samples. Default is 0.
#' @param periodogram_length An integer specifying the number of hours for the Lomb-Scargle periodogram to be performed. Default is NULL and set to the length of \code{time}.
#'
#' @details This function allows for the filtering of luminescence data based on time points and a minimum luminescence threshold in order for the further analysis. Using lomb-scargle periodogram to determine if each sample contains rhythmic data the non-rhythmic data is removed. See X function

#'
#' @return A matrix array is returned with filtered columns and rows based on selected parameters.
#' install.packages("devtools") and install.packages(roxygen2)
#' @import lomb
#' @import stats
#' @import pracma
#' @export
#' @examples
#' data(mydata_example)
#' data <- mydata_example[,-1]
#' time <- as.matrix(mydata_example[,1])
#' o <- filter_curves(data, time, from = 24, to = 144, min_lum = 200, periodogram_length = 48)



filter_curves <- function(x, time, from = min(time, na.rm = TRUE), to = max(time, na.rm = TRUE), min_lum = 0, periodogram_length = NULL) {

  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric. This is not permitted")

  if (!is.matrix(time))
    stop("Time is NOT matrix array. This is not permitted")

  length_df <- length(x)
  rownumber <- function(time, from, to ) {
    minrow <- which(abs(time - from) == min(abs(time - from)))
    maxrow <- which(abs(time - to) == min(abs(time - to)))
    return(c(minrow,maxrow))
  }

  #create data frame for curves and time vector with number of rows based of to and form
  Curves_Timecut <- data.frame(x[rownumber(time = time, from = from, to = to)[1]:rownumber(time = time, from = from, to = to)[2],1:length_df])

  time_cut <- time[rownumber(time = time,from = from,to = to)[1]:rownumber(time = time, from = from, to = to)[2],]

  # create df of curves that exceed the minimum luminescence threshold
  Curves_set_min <- Curves_Timecut[sapply(Curves_Timecut, function(x) min(x, na.rm = T) > min_lum)]

  #trend and then detrend the curves and keep column names the same
  tred <- apply(Curves_set_min, MARGIN = 2, FUN = function(X) stats::lm(X~time_cut, na.action = na.exclude))
  detrend <- as.data.frame(sapply(X = seq(1:length(Curves_set_min)), FUN = function(g){ residuals(tred[[g]])}))
  colnames(detrend) <- colnames(Curves_set_min)

  #depending on length specified to perfom periodogram

  # Calculate length of the periodogram
  if (is.null(periodogram_length)) {
    periodogram_from <- from
    periodogram_to <- to
  } else {
    if (periodogram_length > to - from) {
      stop("Periodogram length exceeds the length between 'from' and 'to'.")
    }
    if (periodogram_length < 18) {
      stop("Periodogram length must be greater than 18")
    }
    periodogram_from <- max(time, na.rm = TRUE) - periodogram_length
    periodogram_to <- max(time, na.rm = TRUE)
  }


  #create dataframe of time from the last number of hours specified by periodogram length
  time_periodogram <- time_cut[rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[1]:rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[2]]


  #create dataframe of the last hours of the curves for periodogram analysis specified by periodogram length
  detrend_periodogram_df <- data.frame(detrend[rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[1]:rownumber(time_cut,max(time_cut) - periodogram_length,max(time_cut))[2],])

  #perform ls periodogram then find the peaks of these and remove any curves where this is below 0.001 threshold
  pero <- apply(detrend_periodogram_df, MARGIN = 2, FUN = function(x) lomb::lsp(x, times = time_periodogram, type =  "period",from = 18, to = 30,ofac = 44, plot = FALSE ))

  peaks <- lapply(X = seq(1:length(Curves_set_min)), FUN = function(g) {pracma::findpeaks(pero[[g]]$power, nups = 1, ndowns = 1, zero = "0")})

  peaks <- lapply(X = seq(1:length(Curves_set_min)), FUN = function(g) {max(peaks[[g]][,1])})

  #this removes curves that do not meet the periodogram threshold
  Curves_filterd <- detrend[sapply(X = seq(1:length(detrend)), FUN = function(g) {peaks[[g]] > pero[[g]]$sig.level})]
  finaldf<-cbind(time_cut,Curves_filterd)
  return(finaldf)

}





