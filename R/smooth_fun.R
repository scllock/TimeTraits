#' Smooth circadian time-series data with optional shape-based outlier detection
#'
#' \code{smooth_fun} smooths circadian time-series data using a functional data
#' analysis (FDA) framework based on cubic B-spline basis functions. Each curve
#' is optionally normalised, smoothed individually, and evaluated on a common
#' time grid. An optional shape-based outlier detection step can be applied using
#' total variation depth, after which the remaining curves are smoothed
#' collectively.
#'
#' The function supports metadata tracking, allowing identification of samples
#' removed during outlier detection.
#'
#' @param x A numeric matrix of observations, where rows correspond to time
#'   points and columns correspond to samples.
#' @param time A numeric matrix giving the time values corresponding to the rows
#'   of \code{x}.
#' @param meta An optional data frame containing sample metadata. Must include a
#'   \code{Sample_id} column matching the column names of \code{x}.
#' @param normalize Logical indicating whether each curve should be normalised
#'   to the range [0, 1] prior to smoothing. Default is \code{TRUE}.
#' @param outlier Logical indicating whether shape-based outlier detection should
#'   be performed using total variation depth. If \code{FALSE}, all curves are
#'   retained and no outlier detection is applied. Default is \code{FALSE}.
#' @param deriv Integer specifying the derivative order at which curves are
#'   evaluated for outlier detection. Default is 1 (first derivative).
#' @param nbasis Integer specifying the maximum number of B-spline basis
#'   functions used when smoothing individual curves. The effective number of
#'   basis functions is constrained by the number of observations per curve.
#'   Default is 50.
#' @param lambda Numeric value controlling the roughness penalty applied during
#'   smoothing of individual curves. Higher values produce smoother curves.
#'   Default is 1.
#'
#' @details
#' Each curve is first optionally normalised and stripped of missing values.
#' Curves are smoothed individually using cubic B-spline basis functions with a
#' second-derivative roughness penalty. Functional parameter objects are cached
#' and reused for curves sharing the same basis configuration to improve
#' computational efficiency.
#'
#' When \code{outlier = TRUE}, curves are smoothed again using a higher roughness
#' penalty and evaluated at the specified derivative order. Shape outliers are
#' identified using total variation depth via
#' \code{fdaoutlier::tvdmss} and removed prior to final smoothing.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{curves}{A numeric matrix of smoothed curves evaluated on a common time
#'     grid, with columns corresponding to retained samples.}
#'   \item{smooth}{A functional data object containing the final smoothed curves.}
#'   \item{time}{A numeric vector giving the common time grid used for evaluation.}
#'   \item{meta_kept}{A data frame of metadata for retained samples, or
#'     \code{NULL} if no metadata were supplied.}
#'   \item{meta_removed}{A data frame of metadata for samples removed during
#'     outlier detection, or \code{NULL} if no samples were removed.}
#' }
#'
#' @import fda
#' @import stats
#' @import fdaoutlier
#' @export
#'
#' #' @examples
#' data(mydata_example)
#' data <- mydata_example[, -1]
#' time <- as.matrix(mydata_example[, 1])
#'
#' smooth_curves <- smooth_fun(
#'   x = data,
#'   time = time,
#'   nbasis = 40,
#'   lambda = 1,
#'   shape = FALSE
#' )
#'
#'
#'
smooth_fun <- function(
    x,
    time,
    meta = NULL,
    normalize = TRUE,
    outlier = FALSE,
    deriv = 1,
    nbasis = 50,
    lambda = 1
)
{

  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric. This is not permitted")

  if (!is.matrix(time))
    stop("Time is NOT matrix array. This is not permitted")

  ## ---- META CHECKS
  if (!is.null(meta)) {
    if (!"Sample_id" %in% colnames(meta))
      stop("meta must contain a Sample_id column")

    if (is.null(colnames(x)))
      stop("x must have column names matching meta$Sample_id")

    if (!all(colnames(x) %in% meta$Sample_id))
      stop("All curve columns must be present in meta$Sample_id")

    meta <- meta[match(colnames(x), meta$Sample_id), ]
  }

  sample_id <- colnames(x)

  ## ---- OPTIONAL NORMALISATION
  if (normalize) {
    x_use <- apply(
      x,
      2,
      function(X)
        (X - min(X, na.rm = TRUE)) / diff(range(X, na.rm = TRUE))
    )
  } else {
    x_use <- x
  }

  ## ---- REMOVE NAs PER CURVE
  everycurve <- vector("list", ncol(x_use))
  for (i in seq_len(ncol(x_use))) {
    combine <- cbind(time, x_use[, i])
    everycurve[[i]] <- na.omit(combine)
  }
  ## ---- BASIS + FDPAR PER CURVE (CACHED)
  sp_fdobj <- vector("list", ncol(x_use))
  fdpar_cache <- new.env(parent = emptyenv())

  for (i in seq_len(ncol(x_use))) {

    curve <- everycurve[[i]]
    nobs  <- nrow(curve)

    if (nobs < 5)
      stop("Each curve must have at least 5 observations for smoothing.")

    nbasis_i <- min(nbasis, nobs - 1)
    key <- paste(nbasis_i, round(min(curve[,1]), 3), round(max(curve[,1]), 3))

    if (!exists(key, fdpar_cache)) {
      basis <- fda::create.bspline.basis(
        rangeval = range(curve[, 1]),
        norder   = 4,
        nbasis   = nbasis_i
      )
      fdpar_cache[[key]] <- fda::fdPar(basis, Lfdobj = 2, lambda = lambda)
    }

    sp_fdobj[[i]] <- fdpar_cache[[key]]
  }

  ## ---- SMOOTH EACH CURVE
  sp_smoothbasis <- vector("list", ncol(x_use))
  for (i in seq_len(ncol(x_use))) {
    sp_smoothbasis[[i]] <- fda::smooth.basis(
      argvals = everycurve[[i]][, 1],
      y       = everycurve[[i]][, 2],
      fdParobj = sp_fdobj[[i]]
    )
  }

  ## ---- COMMON TIME GRID
  min_time <- max(sapply(everycurve, function(x) min(x[, 1])))
  max_time <- min(sapply(everycurve, function(x) max(x[, 1])))
  new_time <- seq(min_time, max_time, length.out = 300)

  ## ---- EVALUATE CURVES
  eval_curves <- vapply(
    sp_smoothbasis,
    function(x) fda::eval.fd(new_time, x$fd),
    FUN.VALUE = numeric(length(new_time))
  )
  colnames(eval_curves) <- sample_id

  Total_curves <- data.matrix(eval_curves)

  ## ---- OUTLIER SMOOTH
  nbasis_outlier <- min(100, length(new_time) - 1)

  BASIS_bspline <- fda::create.bspline.basis(
    rangeval = range(new_time),
    norder   = 4,
    nbasis   = nbasis_outlier
  )


  fdobj_outlier <- fda::fdPar(BASIS_bspline, Lfdobj = 2, lambda = 50)

  sp_totalsmooth_outlier <- fda::smooth.basis(
    argvals = new_time,
    y       = Total_curves,
    fdParobj = fdobj_outlier
  )

  v_sp_totalsmooth_outlier <- fda::eval.fd(
    new_time,
    sp_totalsmooth_outlier$fd,
    deriv
  )


  ## ---- SHAPE OUTLIER DETECTION
  removed_ids <- NULL
  meta_removed <- NULL

  if (outlier) {
    sp_transposed_curves <- t(v_sp_totalsmooth_outlier)

    tvoutlier <- fdaoutlier::tvdmss(dts = sp_transposed_curves)

    if (!is.null(tvoutlier$shape_outliers)) {

      removed_idx <- tvoutlier$shape_outliers

      removed_ids <- sample_id[removed_idx]
      sample_id   <- sample_id[-removed_idx]

      Total_curves <- Total_curves[, -removed_idx, drop = FALSE]

      if (!is.null(meta)) {
        meta_removed <- meta[meta$Sample_id %in% removed_ids, ]
        meta <- meta[meta$Sample_id %in% sample_id, ]
      }
    }
  }

  ## ---- FINAL SMOOTH
  fdobj <- fda::fdPar(BASIS_bspline, Lfdobj = 2, lambda = 1)

  sp_totalsmooth <- fda::smooth.basis(
    argvals = new_time,
    y       = Total_curves,
    fdParobj = fdobj
  )

  colnames(Total_curves) <- sample_id

  ## ---- RETURN
  return(list(
    curves        = Total_curves,
    smooth        = sp_totalsmooth,
    time          = new_time,
    meta_kept     = meta,
    meta_removed  = meta_removed
  ))
}
