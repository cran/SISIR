# Methods for SFCB-class ####
#' @title Methods for SFCB objects
#' @name SFCB-class
#' @export
#' @aliases summary.SFCB
#' @aliases print.SFCB
#' @aliases plot.SFCB
#' @aliases extract_at.SFCB
#' @aliases quality.SFCB
#' @aliases extract_at
#' @aliases quality
#' @aliases SFCB-class
#' @description Print, plot, manipulate or compute quality for outputs of the
#' \code{\link{sfcb}} function (\code{SFCB} object)
#' @param object a \code{SFCB} object
#' @param x a \code{SFCB} object
#' @param ... not used
#' @param plot.type type of the plot. Default to \code{"dendrogram"} (see
#' Details)
#' @param sel.type when \code{plot.type == "selection"}, criterion on which to
#' base the selection. Default to \code{"importance"}
#' @param threshold numeric. When \code{plot.type == "importance"}, numeric
#' threshold to perform a selection (if none has been done before for \code{x})
#' based on thresholding of importance. Default to \code{"none"} in which case
#' no thresholding is performed
#' @param shape.imp when \code{plot.type == "importance"}, type of plot to
#' represent the importance. Default to \code{"boxplot"}
#' @param quality.crit character vector (length 1 or 2) indicating one or two
#' quality criteria to display. The values have to be taken in \{\code{"mse"},
#' \code{"Precision"}, \code{"Recall"}, \code{"ARI"}, \code{"NMI"}\}
#' @param ... not used
#' @param at numeric vector. Set of the number of intervals to extract for 
#' @param ground_truth numeric vector of ground truth. Target variables 
#' to compute qualities correspond to non-zero entries of this vector
#' @param threshold numeric value. If not \code{NULL}, selection of variables to
#' compute qualities is based on a threshold of importance values
#' \code{extract_at}
#' @author {Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' @references {Servien, R. and Vialaneix, N. (2023) A random forest approach 
#' for interval selection in functional regression. Preprint.}
#' @details The \code{plot} functions can be used in four different ways to 
#' extract information from the \code{SFCB} object: \itemize{
#'  \item \code{plot.type == "dendrogram"} displays the dendrogram obtained at
#'  the clustering step of the method. Depending on the cases, the dendrogram
#'  comes with additional information on clusters, variable selections and/or
#'  importance values;
#'  \item \code{plot.type == "selection"} displays either the evolution of the
#'  importance for the simulation with the best (smallest) MSE for each time
#'  step in the range of the functional predictor or the evolution of the
#'  selected intervals along the whole range of the functional prediction also 
#'  for the best MSE;
#'  \item \code{plot.type == "importance"} displays a summary of the importance
#'  values over the whole range of the functional predictor and for the 
#'  different experiments. This summary can take the form of a boxplot or of
#'  an histogram;
#'  \item \code{plot.type == "quality"} displays one or two quality distribution
#'  with respect to the different experiments and different number of intervals.
#' }
#' @seealso \code{\link{sfcb}}
#' @examples 
#' data(truffles)
#' out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
#'              summary.method = "pls", selection.method = "relief")
#' summary(out1)
#' 
#' \dontrun{
#' plot(out1)
#' plot(out1, plot.type = "selection")
#' plot(out1, plot.type = "importance")
#' }
#' 
#' out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
#'              summary.method = "basics", selection.method = "none",
#'              range.at = c(5, 7))
#' out3 <- extract_at(out2, at = 6)
#' summary(out3)
#' 
summary.SFCB <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nSFCB object with:\n")
  lobj <- length(object$groups)
  if (lobj == 1) {
    cat("    -", length(unique(object$groups[[1]])), "interval(s)\n")
  } else {
    minl <- length(unique(object$groups[[1]]))
    maxl <- length(unique(object$groups[[lobj]]))
    cat("    -", minl, "-", maxl, "interval(s)\n")
  }
  if ("selected" %in% names(object)) {
    if (!is.null(ncol(object$selected[[1]]))) {
      selected <- lapply(object$selected, function(alist) {
        apply(alist[, 2:3], 1, sum)
      })
    } else selected <- object$selected
    if (lobj == 1) {
      selected <- object$groups[[1]][selected[[1]]]
      cat("    -", length(unique(selected)), "selected interval(s)\n")
    } else {
      selected <- mapply(function(a, b) a[b], object$groups, selected)
      selected <- sapply(selected, function(x) length(unique(x)))
      cat("    -", min(selected), "-", max(selected), "selected interval(s)\n")
    }
  }
  cat("    -", nrow(object$mse) / lobj, "repeats\n")
  cat("    - MSE ranging in [", min(object$mse$mse), ", ", max(object$mse$mse),
      "]\n", sep = "")
  if ("computational.times" %in% names(object)) {
    cat("    - computational time (total):", sum(object$computational.times), 
        "(seconds)\n\n")
  }
  if ("quality" %in% names(object)) {
    cat("    - precision wrt ground truth in [", 
        min(object$quality$Precision, na.rm = TRUE), ", ", 
        max(object$quality$Precision, na.rm = TRUE), "]\n", sep = "")
    cat("    - recall wrt ground truth in [", 
        min(object$quality$Recall, na.rm = TRUE), ", ", 
        max(object$quality$Recall, na.rm = TRUE), "]\n", sep = "")
    if ("threshold" %in% names(object))
      cat("    for threshold: ", object$threshold, "\n", sep = "")
    cat("\n")
  }
  return(invisible(NULL))
}

#' @export
#' @rdname SFCB-class
print.SFCB <- function(x, ...) {
  summary(x)
}

#' @rdname SFCB-class
#' @export
plot.SFCB <- function(x, ...,
                      plot.type = c("dendrogram", "selection", "importance",
                                    "quality"),
                      sel.type = c("importance", "selection"),
                      threshold = "none", shape.imp = c("boxplot", "histogram"),
                      quality.crit = "mse") {
  args <- list("x" = x)
  plot.type <- match.arg(plot.type)
  args$"shape.imp" <- match.arg(shape.imp)
  args$"sel.type" <- match.arg(sel.type)
  args$"threshold" <- threshold
  args$"quality.crit" <- quality.crit
  
  plot_function <- sprintf("plot_%s", plot.type)
  plot_function <- eval(as.name(plot_function))
  args <- args[names(formals(plot_function))]
  
  p <- do.call("plot_function", args)
  if (plot.type != "dendrogram") return(p)
  
  return(invisible())
}

#' @rdname SFCB-class
#' @export
extract_at <- function(object, at) {
  UseMethod("extract_at")
}

#' @export
extract_at.SFCB <- function(object, at) {
  
  extract_call <- match.call()
  
  if (!is.numeric(at)) stop("'at' must be a numeric vector")
  orig_at <- names(object$groups)
  at <- sort(at)
  selected <- match(at, orig_at)
  if (anyNA(selected)) {
    stop("'at' must be included in the range of tested groups for 'object'.")
  }
  
  out <- object
  out$groups <- out$groups[selected]
  out$summaries <- out$summaries[selected]
  out$mse <- out$mse[out$mse$clust %in% at, ]
  out$importances <- out$importances[selected]
  if ("computational.times" %in% names(out)) out$"computational.times" <- NULL
  if ("selected" %in% names(out)) out$selected <- out$selected[selected]
  
  out$call <- extract_call
  
  return(out)
}

#' @rdname SFCB-class
#' @export
quality <- function(object, ground_truth, threshold = NULL) {
  UseMethod("quality")
}

#' @export
#' @importFrom aricode ARI NID NMI AMI NVI
quality.SFCB <- function(object, ground_truth, threshold = NULL) {
  
  if (length(ground_truth) != length(object$groups[[1]]))
    stop(paste("'ground_truth' must have a length identical to initial number ",
               "of variables."))

  if (!is.null(threshold) && !is.numeric(threshold) && threshold <= 0) 
    stop("'threshold' must be a positive number or NULL.")
  
  if (is.null(threshold) & !("selected" %in% names(object))) 
    stop(paste("No selected interval in this 'SFCB' object and no 'threshold'",
               "provided. Can not compute quality measures..."))
  
  all_at <- unique(object$mse$clust)
  basics <- ncol(object$importances[[1]]) != (nrow(object$mse) / length(all_at))
  out_obj <- object
  
  ground_truth_f <- sapply(ground_truth != 0, ifelse, yes = 1, no = 0)
  ground_truth_f <- factor(ground_truth_f)
  range_names <- names(object$groups)
  var_names <- names(object$groups[[1]])
  if (!is.null(threshold)) {
    if (basics) {
      object$importances <- lapply(object$importances, function(alist) {
        colmean <- grep("mean", colnames(alist))
        colsd <- grep("sd", colnames(alist))
        impm <- lapply(seq_along(colmean), function(ind) {
          pmax(alist[, colmean[ind]], alist[, colsd[ind]], na.rm = TRUE)
        })
        impm <- Reduce(cbind, impm)
        rownames(impm) <- rownames(alist)
        return(impm)
      })
    }
    nbrep <- ncol(object$importances[[1]])
    out <- lapply(object$importances, function(imp) {
      quals <-  lapply(1:ncol(imp), 
                       function(col) rownames(imp)[imp[, col] > threshold])
      quals <- lapply(quals, compute_qualities, var_names, ground_truth_f)
      quals <- Reduce(rbind, quals)
      return(quals)
    })
    out <- data.frame("clust" = rep(range_names, sapply(out, nrow)),
                      "repeats" = rep(1:nbrep, length(all_at)),
                      Reduce(rbind, out))
    rownames(out) <- NULL
  } else {
    if (basics) {
      object$selected <- lapply(object$selected, function(alist) {
        out <- alist$"in.mean" | alist$"in.sd"
        out <- alist$variable[out]
        return(out)
      })
    }
    out <- lapply(object$selected, compute_qualities, var_names = var_names, 
                  ground_truth_f = ground_truth_f)
    out <- data.frame("clust" = range_names, Reduce(rbind, out))
    rownames(out) <- NULL
  }
  
  out_obj$truth <- ground_truth
  out_obj$quality <- out
  if (!is.null(threshold)) out_obj$threshold <- threshold
  return(out_obj)
}

compute_qualities <- function(asel, var_names, ground_truth_f) {
  sel_f <- sapply(var_names %in% asel, ifelse, yes = 1, no = 0)
  sel_f <- factor(sel_f)
  nbpos <- sum(sel_f == "1" & ground_truth_f == "1")
  tpr <- nbpos / sum(sel_f == "1")
  recall <- nbpos / sum(ground_truth_f == "1")
  ari <- ARI(sel_f, ground_truth_f)
  nmi <- NMI(sel_f, ground_truth_f)
  out <- data.frame("Precision" = tpr, "Recall" = recall, "ARI" = ari, 
                    "NMI" = nmi)
  return(out)
}
