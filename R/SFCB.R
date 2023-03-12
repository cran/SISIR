################################################################################
# SFCB
################################################################################
#' @title sfcb
#' @name sfcb
#' @export
#'
#' @description
#' \code{sfcb} performs interval selection based on random forests
#'
#' @param X input predictors (matrix or data.frame)
#' @param Y target variable (vector whose length is equal to the number of rows
#' in X)
#' @param group.method group method. Default to \code{"adjclust"}
#' @param summary.method summary method. Default to \code{"pls"}
#' @param selection.method selection method. Default to \code{"none"} (no 
#' selection performed)
#' @param at number of groups targeted for output results (integer). Not used
#' when \code{range.at} is not \code{NULL}
#' @param range.at (vector of integer) sequence of the numbers of groups for
#' output results
#' @param seed random seed (integer)
#' @param repeats number of repeats for the final random forest computation
#' @param keep.time keep computational times for each step of the method? 
#' (logical; default to \code{TRUE})
#' @param verbose print messages? (logical; default to \code{TRUE})
#' @param parallel not implemented yet
#' 
#' @return an object of class \code{"SFCB"} with elements: 
#'   \item{dendro}{a dendrogram corresponding to the method chosen in
#'   \code{group.method}}
#'   \item{groups}{a list of length \code{length(range.at)} (or of length 1 if
#'   \code{range.at == NULL}) that contains the clusterings of input variables
#'   for the selected group numbers}
#'   \item{summaries}{a list of the same length than \code{$groups} that 
#'   contains the summarized predictors according to the method chosen in 
#'   \code{summary.methods}}
#'   \item{selected}{a list of the same length than \code{$groups} that contains
#'   the names of the variable selected by \code{selection.method} if it is not
#'   equal to \code{"none"}}
#'   \item{mse}{a data.frame with \code{repeats} \eqn{\times}
#'   \code{length($groups)} rows that contains Mean Squared Errors of the 
#'   \code{repeats} random forests fitted for each number of groups}
#'   \item{importance}{a list of the same length than \code{$groups} that 
#'   contains a data.frame providing variable importances for the variables in
#'   selected groups in \code{repeats} columns (one for each iteration of the
#'   random forest method). When \code{summary.method == "basics"}, importance 
#'   for mean and sd are provided in separated columns, in which case, the 
#'   number of columns is equal to 2\code{repeats}}
#'   \item{computational.times}{a vector with 4 values corresponding to the
#'   computational times of (respectively) the group, summary, selection, and RF
#'   steps. Only if \code{keep.time == TRUE}}
#'   \item{call}{function call}
#' 
#' @author {Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' 
#' @references {Servien, R. and Vialaneix, N. (2023) A random forest approach 
#' for interval selection in functional regression. Preprint.}
#' 
#' @examples 
#' data(truffles)
#' out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
#'              summary.method = "pls", selection.method = "relief")
#' out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
#'              summary.method = "basics", selection.method = "none",
#'              range.at = c(5, 7))

sfcb <- function(X, Y, group.method = c("adjclust", "cclustofvar"), 
                 summary.method = c("pls", "basics", "cclustofvar"), 
                 selection.method = c("none", "boruta", "relief"), 
                 at = round(0.15 * ncol(X)), range.at = NULL,  seed = NULL, 
                 repeats = 5, keep.time = TRUE, verbose = TRUE, 
                 parallel = FALSE) {
  
  # Input checking ####
  group.method <- match.arg(group.method)
  summary.method <- match.arg(summary.method)
  selection.method <- match.arg(selection.method)
  if (is.null(range.at)) {
    ## checking 'at'
    at_length_one <- length(at) == 1
    at_numeric <- is.numeric(at)
    if (!at_length_one || !at_numeric) stop("'at' must be a positive integer")
    
    at_positive <- at > 0
    at_integer <- (at - floor(at)) == 0
    if (!at_positive || !at_integer) stop("'at' must be a positive integer")
  } else {
    ## checking 'range.at'
    range_length_two <- length(range.at) == 2
    range_numeric <- is.numeric(range.at)
    if (!range_length_two || !range_numeric) 
      stop(paste("'range.at' must be a vector of minimum and maximum numbers",
                 "of intervals (length 2) or 'NULL'"))
    
    range_positive <- range.at[1] > 0
    range_ordered <- range.at[1] < range.at[2]
    range_integer <- sapply(range.at, function(aat) (aat - floor(aat)) == 0)
    range_integer <- all(range_integer)
    if (!range_positive || !range_ordered || !range_integer) 
      stop(paste("'range.at' must be a vector of minimum and maximum numbers",
                 "of intervals (length 2) or 'NULL'"))
  }
  ## checking 'seed'
  if (!is.null(seed) && !is.numeric(seed)) stop("'seed' must be numeric!")
  
  sfcb_call <- match.call()
  
  # Step1: group computation ####
  group_ct <- system.time({
    if (verbose) message("Creating groups...")
    group_function <- sprintf("group_%s", group.method)
    group_function <- eval(as.name(group_function))
    all_groups <- group_function(X)
  })[3]
  
  # Step2: compute summaries ####
  summary_ct <- system.time({
    if (verbose) message("Computing summaries...")
    if (is.null(range.at)) {
      selected_groups <- all_groups$groups[as.character(at)]
    } else {
      selected <- as.character(range.at[1]:range.at[2])
      selected_groups <- all_groups$groups[selected]
    }
    summary_function <- sprintf("summary_%s", summary.method)
    summary_function <- eval(as.name(summary_function))
    function_inputs <- names(formals(summary_function))
    if ("target" %in% function_inputs) {
      suppressMessages({
        all_summaries <- summary_function(X, Y, selected_groups)
      })
    } else {
      suppressMessages({all_summaries <- summary_function(X, selected_groups)})
    }
    if (summary.method != "basics") {
      all_summaries <- lapply(all_summaries, function(as) {
        out <- as.data.frame(as)
        names(out) <- paste0("g", 1:ncol(as))
        return(out)
      })
    } else {
      all_summaries <- lapply(all_summaries, function(as) {
        out <- as.data.frame(as)
        names(out) <- gsub("clustmean", "mean\\.g", names(out))
        names(out) <- gsub("clustsd", "sd\\.g", names(out))
        return(out)
      })
    }
  })[3]
  
  # Step 3: selection ####
  if (selection.method != "none") {
    selection_ct <- system.time({
      if (verbose) message("Computing selections...")
      
      selection_function <- sprintf("selection_%s", selection.method)
      selection_function <- eval(as.name(selection_function))
      if ("seed" %in% names(formals(selection_function))) {
        all_selections <- selection_function(all_summaries, Y, seed)
      } else all_selections <- selection_function(all_summaries, Y)
      
      # turn selections into their original variables (not at summary level)
      if (summary.method != "basics") {
        orig_selected <- mapply(summary2original, all_selections, 
                                selected_groups, 
                                MoreArgs = list("varnames" = names(X)),
                                SIMPLIFY = FALSE)
      } else {
        orig_selected <- mapply(summary2originalb, all_selections, 
                                selected_groups, all_summaries, 
                                MoreArgs = list("varnames" = names(X)),
                                SIMPLIFY = FALSE)
      }
    })[3]
  } else {
    all_selections <- "none"
    selection_ct <- 0
  }
  
  # Step 4: prediction ####
  rf_ct <- system.time({
    if (verbose) message("Fitting random forest...")
    all_qualities <- predict_rf(all_summaries, Y, selected_groups,
                                all_selections, repeats, seed,
                                (summary.method == "basics"), names(X))
  })[3]
  
  # Output preparation
  if (selection.method != "none") {
    out <- list("dendro" = all_groups$dendro, "groups" = selected_groups,
                "summaries" = all_summaries, "selected" = orig_selected,
                "mse" = all_qualities$mse, 
                "importances" = all_qualities$importances)
  } else {
    out <- list("dendro" = all_groups$dendro, "groups" = selected_groups,
                "summaries" = all_summaries, "mse" = all_qualities$mse, 
                "importances" = all_qualities$importances)
  }
  if (keep.time) {
    out$"computational.times" <- c(group_ct, summary_ct, selection_ct, rf_ct)
    names(out$"computational.times") <- c("group", "summary", "selection", "RF")
  }
  out$call <- sfcb_call
  
  class(out) <- "SFCB"
  
  return(out)
}
