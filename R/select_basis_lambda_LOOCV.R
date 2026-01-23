#' Select optimal number of basis functions and smoothing parameter
#'
#' This function performs leave-one-timepoint-out cross-validation (LOOCV)
#' on a single functional trajectory to determine the optimal number of
#' B-spline basis functions (`nbasis`) and smoothing parameter (`lambda`)
#' for fitting a functional data object using \code{fda::smooth.basis}.
#'
#' Each time point is left out in turn, the curve is refit using the
#' remaining observations, and prediction error at the omitted time
#' point is recorded. Mean squared prediction error is used as the
#' selection criterion.
#'
#' This function is computationally expensive and intended as a helper
#' routine for selecting smoothing parameters prior to large-scale
#' functional analyses.
#' @details
#' The function searches for the best combination of number of basis functions
#' and roughness penalty by minimizing leave-one-out cross-validation error.
#' The smoothing parameter (\eqn{\lambda}) controls the trade-off between
#' goodness of fit and smoothness; larger values enforce smoother curves.

#' @param data A numeric vector of observed values (one curve).
#' @param time A numeric vector of time points corresponding to `data`.
#' @param nbasis_range A numeric vector specifying the range of basis
#'   function numbers to test.
#' @param lambda_range A numeric vector of lambda values to test
#'   (e.g. \code{10^seq(-6, 2, length = 10)}).
#' @param norder Integer specifying the order of the B-spline basis
#'   (default = 4, cubic).
#' @param penalty_order Integer specifying the derivative order used in
#'   the roughness penalty (default = 2).
#'
#' @return A list with:
#' \item{best_nbasis}{Number of basis functions minimizing LOOCV error.}
#' \item{best_lambda}{Lambda value minimizing LOOCV error.}
#' \item{cv_errors}{Matrix of mean squared LOOCV errors across all
#'   combinations.}
#' \item{nbasis_range}{Tested nbasis values.}
#' \item{lambda_range}{Tested lambda values.}
#'
#' @import fda
#' @import stats
#' @export
#'
#' @examples
#' \donttest{
#' time <- seq(0, 48, length.out = 100)
#' y <- sin(2 * pi * time / 24) + rnorm(100, 0, 0.1)
#' result <- select_basis_lambda_LOOCV(y, time,
#'                                     nbasis_range = seq(5, 25, 5),
#'                                     lambda_range = 10^seq(-6, 0, length = 5))
#' print(result)}
#'
select_basis_lambda_LOOCV <- function(
    data,
    time,
    nbasis_range = seq(10, 100, 5),
    lambda_range = 10^seq(-6, 2, length.out = 10),
    norder = 4,
    penalty_order = 2
) {

  # ------------------ checks ------------------
  if (!is.numeric(data)) stop("data must be numeric.")
  if (!is.numeric(time)) stop("time must be numeric.")
  if (length(data) != length(time))
    stop("data and time must have the same length.")
  if (penalty_order >= norder)
    stop("penalty_order must be less than norder.")

  # enforce ordering
  ord <- order(time)
  time <- time[ord]
  data <- data[ord]

  m <- length(time)

  if (m > 200)
    warning("LOOCV is computationally expensive; consider subsampling time points.")

  cv_errors <- matrix(
    NA_real_,
    nrow = length(nbasis_range),
    ncol = length(lambda_range),
    dimnames = list(
      paste0("nbasis_", nbasis_range),
      paste0("lambda_", lambda_range)
    )
  )

  # ------------------ LOOCV grid search ------------------
  for (i in seq_along(nbasis_range)) {

    nbasis <- nbasis_range[i]

    for (j in seq_along(lambda_range)) {

      lambda <- lambda_range[j]
      SSE <- rep(NA_real_, m)

      for (k in seq_len(m)) {

        time_train <- time[-k]
        y_train    <- data[-k]

        # create basis on training range (important for stability)
        basis <- create.bspline.basis(
          rangeval = range(time),
          nbasis   = nbasis,
          norder   = norder
        )

        fdParobj <- fdPar(
          basis,
          Lfdobj = int2Lfd(penalty_order),
          lambda = lambda
        )

        fit <- try(
          smooth.basis(time_train, y_train, fdParobj),
          silent = TRUE
        )

        if (inherits(fit, "try-error"))
          next

        y_pred <- eval.fd(time[k], fit$fd)
        SSE[k] <- (data[k] - y_pred)^2
      }

      cv_errors[i, j] <- mean(SSE, na.rm = TRUE)
    }
  }

  # ------------------ extract optimum ------------------
  best_idx <- which(cv_errors == min(cv_errors, na.rm = TRUE), arr.ind = TRUE)

  best_nbasis <- nbasis_range[best_idx[1, 1]]
  best_lambda <- lambda_range[best_idx[1, 2]]

  list(
    best_nbasis  = best_nbasis,
    best_lambda  = best_lambda,
    cv_errors    = cv_errors,
    nbasis_range = nbasis_range,
    lambda_range = lambda_range
  )
}
