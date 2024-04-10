#' Compute and filter circadian time series data
#'
#' The function is designed to smooth multiple curves represented by time-series data. It normalises each curve and organizes them into a list, removing any missing values. It then creates basis functions and functional parameters for each curve using B-splines the default is set to the number of timepoints and sets a roughness penalty. Next, it generates smooth basis functions for each curve individually. After defining a new time vector covering the range of all curves, it evaluates the smoothed curves under this new time vector. Subsequently, it combines all evaluated curves and performs an outlier detection process using the Total Variation (TV) smooth method from \code{fdaoutliers}. Any shape outliers are removed, and the remaining curves are smoothed again collectively. Finally, the function returns the smoothed curves along with the new time vector and the total curves matrix.
#'
#'
#' @param x Matrix containing data values, each row is a time point and each column corresponds to a sample.
#' @param time matrix containing the time values.
#' @param shape A logical indicating whether shape outlier detection should be performed
#' @param outlier A positive integer specifying the derivative the curves will be evaluated at for outlier detection.The default is the first derivative.
#'
#' @details This function allows estimation of functions of time series data using cubic b-spline basis system.

#'
#' @return A smooth functional data object
#'
#' @import fda
#' @import stats
#' @import fdaoutlier
#' @export
#' @examples
#' data(mydata_example)
#' data <- mydata_example[,-1]
#' time <- as.matrix(mydata_example[,1])
#' smooth_curves <- smooth_fun(data, time, shape = FALSE, outlier = 1)
#'
#'
#'

smooth_fun <- function(x, time, outlier = TRUE, shape = TRUE) {

  if (!any(apply(x, 2, is.numeric)))
    stop("values are NOT numeric. This is not permitted")

  if (!is.matrix(time))
    stop("Time is NOT matrix array. This is not permitted")
  # NORMALISE each curve
  norm_all <- apply(x, MARGIN = 2, FUN = function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
  #create list containing each column and the time column with NAs removed
  everycurve <- list()
  for (i in 1:ncol(norm_all)) {

    #combine time column with every every column individually
    combine <- cbind(time,norm_all[,i])
    #append list with each time and curve combo with NAs removed.
    everycurve[[i]] <- na.omit(combine)
  }
  #now make a basis functions and functional parameter for each one....
  sp_fdobj <- list()
  for (i in 1:ncol(norm_all)) {
    #range and number of basis for each bassis functtions
    range <- c(min(everycurve[[i]][,1]), max(everycurve[[i]][,1]))
    basis_num <- length(everycurve[[i]][,1]) + 2
    # basis function
    sp_BASIS_bspline <- fda::create.bspline.basis(rangeval = range ,norder = 4,nbasis = basis_num)
    #Define Functional Parameter with Roughness Penalty
    sp_fdobj[[i]] <- fda::fdPar(sp_BASIS_bspline,Lfdobj = 2,lambda = 1)
  }
  # for each curve make the smoothbasis
  sp_smoothbasis <- list()
  for (i in 1:ncol(norm_all)) {
    sp_smoothbasis[[i]] <- fda::smooth.basis(argvals = everycurve[[i]][,1], y = everycurve[[i]][,2], fdParobj = sp_fdobj[[i]])
  }
  #make a new time vector that will cover all the different times
  #the max of the minamum ZT time and the minamum of the maximum time
  min_time <- max(sapply(everycurve, function(x) min(x[,1])))
  max_time <- min(sapply(everycurve, function(x) max(x[,1])))
  new_time <- seq(min_time,max_time,length.out = 300)

  #evaluate the curves all under this new time
  eval_curves <- sapply(sp_smoothbasis, function(x) fda::eval.fd(new_time,x$fd))
  colnames(eval_curves) <- colnames(norm_all)

  #now create smooth functions for all cuvres together!
  Total_curves <- data.matrix(eval_curves)
  BASIS_bspline <- fda::create.bspline.basis(rangeval = c(min(new_time),max(new_time)),norder = 4,nbasis = 302)
  fdobj_outlier = fda::fdPar(BASIS_bspline,Lfdobj = 2,lambda = 50)
  sp_totalsmooth_outlier <- fda::smooth.basis(argvals = new_time, y = Total_curves, fdParobj = fdobj_outlier)


  v_sp_totalsmooth_outlier <- fda::eval.fd(new_time,sp_totalsmooth_outlier$fd, outlier)
  #transpose and do tvs mss outlier detection....

  if (shape) {
    sp_transposed_curves <- data.matrix(t(v_sp_totalsmooth_outlier))


    tvoutlier <- fdaoutlier::tvdmss(dts = sp_transposed_curves)

    #tvoutlier$shape_outliers
    #tvoutlier$magnitude_outliers
    if (!is.null(tvoutlier$shape_outliers)) {
      Total_curves <- data.matrix(Total_curves[, -tvoutlier$shape_outliers])
    }
  }
  fdobj = fda::fdPar(BASIS_bspline,Lfdobj = 2,lambda = 1)
  sp_totalsmooth <- fda::smooth.basis(argvals = new_time, y = Total_curves, fdParobj = fdobj)

  return(sp_totalsmooth)
  return(new_time)
  return(Total_curves)
}





# create basis function
# CV estimation
# create functional data object
# FDA outlier.....









