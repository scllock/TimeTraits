#' Filter, detrend, and quality-control luminescence time-series data
#'
#' \code{filter_curves} preprocesses luminescence time-series data by applying
#' a series of quality control and rhythm-detection steps. The function first
#' subsets the data to a specified time window, removes samples with low signal
#' intensity, detrends each time series using linear regression, and then
#' assesses rhythmicity using Lomb–Scargle periodogram analysis. Only samples
#' showing statistically significant rhythmic components are retained.
#'
#' The function returns detrended curves that pass all filtering criteria,
#' along with optional metadata for retained samples and a log of samples
#' removed at each step.
#'
#' @param x A numeric matrix of luminescence values, where rows correspond to
#'   time points and columns correspond to samples. Sparse matrices are allowed.
#' @param time A numeric vector giving the time points corresponding to the rows
#'   of \code{x}.
#' @param meta An optional data frame containing sample metadata. Must include
#'   a \code{Sample_id} column and may include a \code{Genotype} column.
#' @param from Numeric value specifying the start time for analysis. Defaults
#'   to the minimum of \code{time}.
#' @param to Numeric value specifying the end time for analysis. Defaults to
#'   the maximum of \code{time}.
#' @param min_lum Numeric threshold specifying the minimum luminescence value
#'   required for a sample to be retained. Samples with a minimum value less
#'   than or equal to this threshold are removed. Default is 0.
#' @param periodogram_length Numeric value specifying the length (in the same
#'   units as \code{time}) of the time window used for Lomb–Scargle periodogram
#'   analysis. If \code{NULL}, the full selected time range is used.
#'
#' @details
#' After subsetting the data to the requested time window, samples are filtered
#' based on minimum luminescence to remove low-signal or failed measurements.
#' Each remaining time series is detrended by fitting a linear model against
#' time and extracting residuals. Rhythmicity is then assessed using a
#' Lomb–Scargle periodogram over periods between 18 and 30 hours. Samples whose
#' strongest spectral peak does not exceed the corresponding significance
#' threshold are classified as non-rhythmic and removed.
#'
#' A record of all removed samples and the reason for their removal is returned
#' to support transparent quality control.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{curves}{A data frame containing the time vector and detrended
#'     luminescence curves that passed all filtering steps.}
#'   \item{meta_kept}{A data frame of metadata corresponding to retained samples,
#'     or \code{NULL} if no metadata were supplied.}
#'   \item{removed}{A data frame logging samples removed during filtering and
#'     the reason for removal, or \code{NULL} if no samples were removed.}
#' }
#'
#' @import stats
#' @import lomb
#' @import pracma
#' @export
#'
#' @examples
#' data(mydata_example)
#' data <- mydata_example[, -1]
#' time <- as.numeric(mydata_example[, 1])
#'
#' res <- filter_curves(
#'   x = data,
#'   time = time,
#'   from = 24,
#'   to = 144,
#'   min_lum = 200,
#'   periodogram_length = 48
#' )
filter_curves <- function(
    x,
    time,
    meta = NULL,
    from = min(time, na.rm = TRUE),
    to   = max(time, na.rm = TRUE),
    min_lum = 0,
    periodogram_length = NULL
) {
  
  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric")
  
  if (!is.numeric(time))
    stop("Time must be a numeric vector")
  
  ## ---- helpers ----
  rownumber <- function(time, from, to) {
    c(
      which.min(abs(time - from)),
      which.min(abs(time - to))
    )
  }
  
  idx <- rownumber(time, from, to)
  
  Curves_Timecut <- as.data.frame(x[idx[1]:idx[2], , drop = FALSE])
  time_cut <- time[idx[1]:idx[2]]
  
  sample_ids <- colnames(Curves_Timecut)
  if (is.null(sample_ids))
    sample_ids <- paste0("curve_", seq_len(ncol(Curves_Timecut)))
  
  ## ---- initialise removal log ----
  removed_df <- data.frame(
    Sample_id = character(),
    Genotype  = character(),
    Reason    = character(),
    stringsAsFactors = FALSE
  )
  
  ## ---- min lum filter ----
  min_lum_pass <- sapply(Curves_Timecut, function(z)
    min(z, na.rm = TRUE) > min_lum
  )
  
  if (any(!min_lum_pass)) {
    removed_df <- rbind(
      removed_df,
      data.frame(
        Sample_id = sample_ids[!min_lum_pass],
        Genotype  = NA,
        Reason    = "min_lum",
        stringsAsFactors = FALSE
      )
    )
  }
  
  Curves_set_min <- Curves_Timecut[, min_lum_pass, drop = FALSE]
  if (ncol(Curves_set_min) == 0)
    return(list(curves = NULL, removed = removed_df))
  
  ## ---- detrend ----
  tred <- apply(Curves_set_min, 2, function(y)
    lm(y ~ time_cut, na.action = na.exclude)
  )
  
  detrend <- as.data.frame(
    sapply(seq_along(tred), function(i) residuals(tred[[i]]))
  )
  colnames(detrend) <- colnames(Curves_set_min)
  
  ## ---- periodogram window ----
  if (is.null(periodogram_length))
    periodogram_length <- length(time_cut)
  
  idxp <- rownumber(
    time_cut,
    max(time_cut) - periodogram_length,
    max(time_cut)
  )
  
  time_p <- time_cut[idxp[1]:idxp[2]]
  det_p  <- detrend[idxp[1]:idxp[2], , drop = FALSE]
  
  ## ---- lomb–scargle ----
  pero <- apply(det_p, 2, function(y)
    lomb::lsp(
      y,
      times = time_p,
      type = "period",
      from = 18,
      to   = 30,
      ofac = 44,
      plot = FALSE
    )
  )
  
  peaks <- sapply(seq_along(pero), function(i) {
    pk <- pracma::findpeaks(pero[[i]]$power, nups = 1, ndowns = 1)
    if (is.null(pk)) return(NA_real_)
    max(pk[, 1])
  })
  
  keep_idx <- !is.na(peaks) &
    peaks > sapply(pero, `[[`, "sig.level")
  
  if (any(!keep_idx)) {
    removed_df <- rbind(
      removed_df,
      data.frame(
        Sample_id = colnames(detrend)[!keep_idx],
        Genotype  = NA,
        Reason    = "periodogram",
        stringsAsFactors = FALSE
      )
    )
  }
  
  Curves_filtered <- detrend[, keep_idx, drop = FALSE]
  if (ncol(Curves_filtered) == 0)
    return(list(curves = NULL, removed = removed_df))
  
  ## ---- attach metadata if present ----
  if (!is.null(meta)) {
    removed_df <- merge(
      removed_df,
      meta[, c("Sample_id", "Genotype")],
      by = "Sample_id",
      all.x = TRUE,
      suffixes = c("", ".meta")
    )
    removed_df$Genotype <- ifelse(
      is.na(removed_df$Genotype.meta),
      removed_df$Genotype,
      removed_df$Genotype.meta
    )
    removed_df$Genotype.meta <- NULL
  }
  #meta data 
  meta_kept <- NULL
  if (!is.null(meta)) {
    kept_ids <- colnames(Curves_filtered)
    meta_kept <- meta[meta$Sample_id %in% kept_ids, , drop = FALSE]
    meta_kept <- meta_kept[match(kept_ids, meta_kept$Sample_id), ]
  }
  
  ## ---- final output ----
  finaldf <- cbind(time = time_cut, Curves_filtered)
  
  list(
    curves    = finaldf,
    meta_kept = meta_kept,
    removed   = if (nrow(removed_df) > 0) removed_df else NULL
  )
}



