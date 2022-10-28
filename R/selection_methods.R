#' @importFrom Boruta Boruta
#' @importFrom CORElearn attrEval

selection_boruta <- function(summaries, Y, seed) {
  # Description #####
  # returns a list with selected indices for each level of the hierarchy of
  # groups and summaries (Boruta method)
  # Inputs #####
  # summaries: a list with summaries obtained for each level of the hierarchy of
  # groups
  # Y: corresponding target
  # Outputs #####
  # a list with selected indices (same length as summaries)
  
  selections <- sapply(summaries, simplify = FALSE, function(curX) {
    all_indices <- 1:ncol(curX)
    has_nona <- colSums(is.na(curX)) == 0
    curX <- curX[, has_nona, drop = FALSE]
    cur_df <- data.frame(curX, "Y" = Y)
    if (!is.null(seed)) set.seed(seed)
    boruta_decision <- Boruta::Boruta(Y ~ ., data = cur_df, num.threads = 1)
    boruta_decision <- boruta_decision$finalDecision
    out <- which(boruta_decision %in% c("Tentative", "Confirmed"))
    out <- all_indices[has_nona][out]
    out <- unname(out)
    
    return(out)
  })
  
  return(selections)
}

selection_relief <- function(summaries, Y) {
  # Description #####
  # returns a list with selected indices for each level of the hierarchy of
  # groups and summaries (Relief method)
  # Inputs #####
  # summaries: a list with summaries obtained for each level of the hierarchy of
  # groups
  # Y: corresponding target
  # Outputs #####
  # a list with selected indices (same length as summaries)
  
  selections <- sapply(summaries, simplify = FALSE, function(curX) {
    all_indices <- 1:ncol(curX)
    has_nona <- colSums(is.na(curX)) == 0
    curX <- curX[, has_nona, drop = FALSE]
    df <- data.frame("y" = Y, curX)
    eval_f <- CORElearn::attrEval(y ~ ., data = df, 
                                  estimator = "RReliefFexpRank")
    evol <- eval_f[order(eval_f, decreasing = TRUE)] + min(eval_f)
    n <- ncol(df) - 1
    nb_pos <- max(c(1, sum(eval_f > 0)))
    tot_evol <- sum(evol)
    bs <- tot_evol * rev(cumsum(1/(n:1))/n)
    cdt <- which(evol <= bs)
    if (length(cdt) > 0) {
      nb <- min(c(min(cdt), nb_pos))
    } else nb <- nb_pos
    ordered_f <- length(eval_f) + 1 - rank(eval_f)
    out <- which(ordered_f <= nb)
    out <- all_indices[has_nona][out]
    out <- unname(out)
    
    return(out)
  })
  
  return(selections)
}