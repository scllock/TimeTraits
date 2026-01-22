#' Functional PCA traits for circadian rhythm data
#'
#' Performs functional principal component analysis (FPCA) on smoothed
#' circadian time-series data and extracts group-level traits from the
#' resulting FPCA score space.
#'
#' The analysis can optionally be stratified by time segments
#' (e.g. pre/post environmental shift) and by derivatives of the curves.
#' When segments are supplied, functional observations from all segments
#' are first re-evaluated onto a shared time grid and analysed jointly,
#' ensuring that FPCA scores are directly comparable across segments.
#'
#' For each derivative, FPCA is performed once across all samples and
#' segments simultaneously. Segment- and group-specific traits are then
#' computed post hoc from the shared FPCA score space.
#'
#' Extracted traits include:
#' \itemize{
#'   \item Mean and standard deviation of FPCA scores (per principal component)
#'   \item Mean and standard deviation of distance to the FPCA-space centroid
#'   \item Convex hull area in PC1--PC2 space (measure of within-group diversity)
#' }
#'
#' @param smooth_out An object of class \code{fdSmooth}, typically the output
#'   of a smoothing function applied to circadian time-series data.
#'
#' @param time A numeric vector giving the time points corresponding to the
#'   functional observations.
#'
#' @param meta An optional data frame containing sample metadata. Must include
#'   a column named \code{Sample_id} and a grouping column specified by
#'   \code{group_col}.
#'
#' @param group_col Character string giving the name of the column in
#'   \code{meta} that defines the grouping variable (e.g. genotype).
#'   Default is \code{"Genotype"}.
#'
#' @param groups Optional factor or vector defining group membership for each
#'   sample. Used only if \code{meta} is not supplied. Length must match the
#'   number of curves.
#'
#' @param derivative Integer vector specifying which derivatives of the
#'   functional data to analyse (e.g. 0 = original curve, 1 = first derivative).
#'   Default is \code{c(0, 1, 2)}.
#'
#' @param n_pc Integer giving the number of principal components to retain.
#'   If \code{NULL}, the maximum number of PCs allowed by the data is used.
#'
#' @param segments Optional named list defining time segments over which
#'   traits should be summarised. Each element must be a numeric vector of
#'   length two giving the start and end times (e.g.
#'   \code{list(pre = c(0, 104), post = c(104, 240))}).
#'   If \code{NULL}, the full time range is used.
#'
#' @param centerfns Logical indicating whether to centre functional observations
#'   before FPCA. Passed to \code{\link[fda]{pca.fd}}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{traits}{A tidy data frame of group-level FPCA-derived traits,
#'     including marginal PC summaries and FPCA-space dispersion measures.}
#'   \item{scores}{A data frame of FPCA scores for each sample, annotated with
#'     group, segment, and derivative information.}
#'   \item{fpca}{A named list of \code{\link[fda]{pca.fd}} objects, one per
#'     derivative.}
#'   \item{meta}{The metadata data frame aligned to the FPCA input order.}
#' }
#'
#' @details
#' FPCA is performed in a shared functional space for each derivative.
#' When segments are provided, all segment-specific curves are interpolated
#' onto a common time grid prior to FPCA, ensuring statistical comparability
#' of FPCA scores across segments. Group-level traits are computed after FPCA
#' and do not influence the decomposition itself.
#'
#' Convex hull area is computed in PC1--PC2 space only and requires at least
#' three samples per group and segment.
#'
#' @seealso \code{\link[fda]{pca.fd}}, \code{\link[fda]{fdSmooth}},
#'   \code{\link[fda]{eval.fd}}
#' @export
#' @examples
#' \dontrun{
#' fpca_res <- functional_traits(
#'   smooth_out = smooth_data,
#'   time       = time_vec,
#'   meta       = meta_df,
#'   n_pc       = 3,
#'   derivative = c(0, 1, 2),
#'   segments   = list(
#'     pre  = c(56, 104),
#'     post = c(104, 152)
#'   )
#' )
#'
#' head(fpca_res$traits)
#' }
#'

functional_traits <- function(
    smooth_out,
    time,
    meta = NULL,
    group_col = "Genotype",
    groups = NULL,
    derivative = c(0, 1, 2),
    n_pc = NULL,
    segments = NULL,
    centerfns = FALSE) 
 {
  
  # ------------------ checks ------------------
  if (!inherits(smooth_out, "fdSmooth"))
    stop("smooth_out must be smooth.basis output from smooth_fun()")
  
  sample_ids <- colnames(smooth_out$fd$coefs)
  n_curves <- length(sample_ids)
  
  # ------------------ metadata ------------------
  if (!is.null(meta)) {
    
    if (!all(c("Sample_id", group_col) %in% colnames(meta)))
      stop("meta must contain Sample_id and the specified group_col")
    
    meta <- meta[match(sample_ids, meta$Sample_id), ]
    
    if (any(is.na(meta$Sample_id)))
      stop("Sample_id mismatch between smooth_out and meta")
    
    groups <- as.factor(meta[[group_col]])
    
  } else {
    
    if (is.null(groups))
      stop("Either meta or groups must be provided")
    
    if (length(groups) != n_curves)
      stop("Length of groups must match number of curves")
    
    groups <- as.factor(groups)
    meta <- data.frame(Sample_id = sample_ids, group = groups)
  }
  
  # ------------------ segments ------------------
  if (is.null(segments))
    segments <- list(all = range(time))
  
  # ------------------ helpers ------------------
  segment_eval <- function(fd, seg, time) {
    idx <- which(time >= seg[1] & time <= seg[2])
    list(
      time = time[idx],
      values = fda::eval.fd(time[idx], fd)
    )
  }
  
  hull_area <- function(x, y) {
    if (length(x) < 3) return(NA_real_)
    idx <- chull(x, y)
    poly <- cbind(x[idx], y[idx])
    abs(sum(
      poly[-1,1] * poly[-nrow(poly),2] -
        poly[-nrow(poly),1] * poly[-1,2]
    )) / 2
  }
  
  # ================== EVALUATE ALL SEGMENTS ==================
  eval_list <- list()
  
  for (seg_name in names(segments)) {
    eval_list[[seg_name]] <- segment_eval(
      smooth_out$fd,
      segments[[seg_name]],
      time
    )
  }
  
  # ------------------ COMBINE SEGMENTS ------------------
  common_time <- time
  
  all_values <- fda::eval.fd(common_time, smooth_out$fd)
  
  # replicate curves per segment
  all_values <- do.call(cbind, replicate(
    length(segments),
    all_values,
    simplify = FALSE
  ))
  
  all_segments <- rep(names(eval_list), each = n_curves)
  all_sample_ids <- rep(sample_ids, times = length(eval_list))
  all_groups <- rep(groups, times = length(eval_list))
  
  # ------------------ RE-SMOOTH COMBINED DATA ------------------
  basis <- fda::create.bspline.basis(
    rangeval = range(common_time),
    norder = 4,
    nbasis = min(length(common_time) - 1, 50)
  )
  
  fd_all <- fda::smooth.basis(
    argvals = common_time,
    y = all_values,
    fdParobj = fda::fdPar(basis, Lfdobj = 2, lambda = 1)
  )$fd
  # ================== FPCA (ONCE PER DERIVATIVE) ==================
  trait_list <- list()
  score_list <- list()
  fpca_list  <- list()
  
  for (d in derivative) {
    
    fd_deriv <- fda::deriv.fd(fd_all, d)
    
    fpca <- fda::pca.fd(
      fd_deriv,
      nharm = n_pc,
      centerfns = centerfns
    )
    
    fpca_list[[paste0("deriv_", d)]] <- fpca
    
    scores <- as.data.frame(fpca$scores)
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    
    scores$Sample_id <- all_sample_ids
    scores$group <- all_groups
    scores$segment <- all_segments
    scores$derivative <- d
    
    score_list[[paste0("deriv_", d)]] <- scores
    
    # ------------------ PC summary traits ------------------
    pc_traits <- do.call(
      rbind,
      lapply(seq_len(ncol(scores) - 4), function(pc) {
        
        agg <- aggregate(
          scores[[pc]],
          by = list(
            group = scores$group,
            segment = scores$segment,
            derivative = scores$derivative
          ),
          FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x))
        )
        
        do.call(
          rbind,
          lapply(seq_len(nrow(agg)), function(i) {
            data.frame(
              group      = agg$group[i],
              segment    = agg$segment[i],
              derivative = agg$derivative[i],
              trait      = paste0("PC", pc),
              stat       = c("mean", "sd"),
              value      = c(agg$x[i, "mean"], agg$x[i, "sd"]),
              n          = agg$x[i, "n"]
            )
          })
        )
      })
    )
  }
  # ------------------ GEOMETRY TRAITS ------------------
  geom_traits <- do.call(
    rbind,
    split(scores, list(scores$group, scores$segment, scores$derivative),
          drop = TRUE) |>
      lapply(function(df) {
        
        out <- list()
        
        for (pc in paste0("PC", seq_len(ncol(scores) - 4))) {
          mu <- mean(df[[pc]])
          d0 <- abs(df[[pc]] - mu)
          
          out[[pc]] <- data.frame(
            group      = df$group[1],
            segment    = df$segment[1],
            derivative = df$derivative[1],
            trait      = paste0(pc, "_dist_centroid"),
            stat       = c("mean", "sd"),
            value      = c(mean(d0), sd(d0)),
            n          = length(d0)
          )
        }
        
        hull <- if (all(c("PC1", "PC2") %in% colnames(df))) {
          data.frame(
            group      = df$group[1],
            segment    = df$segment[1],
            derivative = df$derivative[1],
            trait      = "PC1_PC2_hull_area",
            stat       = "value",
            value      = hull_area(df$PC1, df$PC2),
            n          = nrow(df)
          )
        }
        
        do.call(rbind, c(out, list(hull)))
      })
  )
  trait_list[[paste0("deriv_", d)]] <- rbind(
    pc_traits,
    geom_traits
  )
  # ================== RETURN ==================
  list(
    traits = do.call(rbind, trait_list),
    scores = do.call(rbind, score_list),
    fpca   = fpca_list,
    meta   = meta
  )
}