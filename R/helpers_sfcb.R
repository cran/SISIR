# related to group selection ####

summary2original <- function(selected, groups, varnames) {
  input_selected <- varnames[groups %in% selected]
  return(input_selected)
}

summary2originalb <- function(selected, groups, summaries, varnames) {
  input_selected <- names(summaries)[selected]
  mean_selected <- grep("mean", input_selected, value = TRUE)
  mean_selected <- as.numeric(gsub("mean\\.g", "", mean_selected))
  sd_selected <- grep("sd", input_selected, value = TRUE)
  sd_selected <- as.numeric(gsub("sd\\.g", "", sd_selected))
  input_selected <- varnames[groups %in% unique(c(mean_selected, sd_selected))]
  in_mean <- varnames[groups %in% mean_selected]
  in_sd <- varnames[groups %in% sd_selected]
  input_selected <- data.frame("variable" = input_selected, 
                               "in.mean"  = input_selected %in% in_mean,
                               "in.sd"    = input_selected %in% in_sd)
  return(input_selected)
}

# related to random forest fitting ####
#' @importFrom ranger ranger

predict_rf <- function(summaries, target, input_groups, selection = "none", 
                       repeats, seed, is_basics = FALSE, varnames) {
  # Description #####
  # returns a list with information on qualities of RF (importance, MSE) for the
  # given summary list and groups to predict a target
  # Inputs #####
  # summaries: a list with a hierarchy of summaries
  # target: the variable to predict
  # input_groups: a list with a hierarchy of groups corresponding to summaries
  # selection (optional): a list with selected features at each level of the
  # hierarchy
  # repeats: number of times the RF is fitted
  # seed: random seed
  # is_basics: are the simulations based on basics predictions
  # varnames: names of input variables
  # Outputs #####
  # a list with two entries: 1/ qualities is a data.frame where rows correspond
  # to experiments (5 for each level of the hierarchy) and columns contain 
  # MSE, number of clusters and (optionally) mean importance for truly 
  # important and non important variables as well as p-value of the Wilcoxon
  # test comparing the distribution of importances in the two groups of 
  # variables; 2/ importance is a list where each entry corresponds to an 
  # experiment and contains importance of the original variables as well as 
  # (optionally) importance of the truly important and non important variables
  
  summary_names <- names(summaries)
  
  if (length(selection) > 1 || selection != "none") {
    summaries <- lapply(1:length(summaries), function(ind) {
      out <- summaries[[ind]][, selection[[ind]], drop = FALSE]
      colnames(out) <- colnames(summaries[[ind]])[selection[[ind]]]
      return(out)
    })
  }
  
  ## compute random forests
  quality_criteria <- lapply(summaries, compute_rfs, target = target, 
                             repeats = repeats, seed = seed)
  quality_criteria <- Reduce(c, quality_criteria)
  
  ## MSE
  all_mse <- sapply(quality_criteria, function(alist) alist$mse)
  all_mse <- data.frame("clust" = rep(summary_names, each = repeats),
                        "mse" = all_mse)
  all_mse$clust <- as.numeric(all_mse$clust)
  
  ## extracting importance
  all_importances <- sapply(quality_criteria, function(alist) alist$importance,
                            simplify = FALSE)
  freq_groups <- sapply(input_groups, table, simplify = FALSE)
  orig_importance <- lapply(1:length(all_importances), compute_origimportance,
                            importances  = all_importances, 
                            input_groups = input_groups, 
                            freq_groups  = freq_groups,
                            varnames     = varnames,
                            repeats      = repeats, 
                            is_basics    = is_basics)
  orig_importance <- sapply(1:length(summaries), function(ind) {
    out <- Reduce(cbind, orig_importance[(repeats*(ind-1)+1):(repeats*ind)])
    if (!is_basics) {
      colnames(out) <- paste0("rep", 1:repeats)
    } else {
      colnames(out) <- paste0(rep(c("mean", "sd"), repeats), 
                              rep(1:repeats, each = 2))
    }
    return(out)
  }, simplify = FALSE)
  names(orig_importance) <- summary_names
  
  out <- list("mse" = all_mse, "importances" = orig_importance)
  
  return(out)
}

compute_rfs <- function(summaries, target, repeats, seed) {
  if (!any(is.na(summaries)) && ncol(summaries) > 0) {
    
    curdf <- data.frame(y = target, summaries)
    
    set.seed(seed)
    random_seeds <- sample(1:10^6, repeats, replace = FALSE)
    cur_model <- sapply(1:repeats, function(ind) {
      ranger::ranger(y ~ ., data = curdf, importance = "permutation", 
                     seed = random_seeds[ind], num.threads = 6)
    }, simplify = FALSE)
    out <- lapply(cur_model, function(ares) {
      list("mse" = ares$prediction.error,
           "importance" = ares$variable.importance)
    })
    
  } else {
    
    out <- lapply(1:repeats, function(ind) 
      list("mse" = NA, "importance" = rep(NA, ncol(summaries))))
  
  }
  
  return(out)
}

compute_origimportance <- function(ind, importances, input_groups, freq_groups,
                                   varnames, repeats, is_basics) {
  expe <- (ind - 1)%/%repeats + 1
    
  if (!is_basics) {
    
    groups_in_imp <- as.numeric(gsub("g", "", names(importances[[ind]])))
    frequencies <- freq_groups[[expe]][as.character(groups_in_imp)]
    out <- rep(importances[[ind]], frequencies)
    names(out) <- varnames[input_groups[[expe]] %in% groups_in_imp]
  
  } else {
    
    sel_in_means <- grep("mean", names(importances[[ind]]), value = TRUE)
    groups_in_means <-  gsub("mean\\.g", "", sel_in_means)
    sel_in_sd <- grep("sd", names(importances[[ind]]), value = TRUE)
    groups_in_sd <-  gsub("sd\\.g", "", sel_in_sd)
    freq_means <- freq_groups[[expe]][groups_in_means]
    freq_sd <- freq_groups[[expe]][groups_in_sd]
    out_means <- rep(importances[[ind]][sel_in_means], freq_means)
    names_means <- varnames[input_groups[[expe]] %in% groups_in_means]
    out_sd <- rep(importances[[ind]][sel_in_sd], freq_sd)
    names_sd <- varnames[input_groups[[expe]] %in% groups_in_sd]
    out <- merge(data.frame("var" = names_means, "mean" = out_means), 
                 data.frame("var" = names_sd, "sd" = out_sd),
                 by = "var", all = TRUE)
    rownames(out) <- out$var
    out <- out[varnames, ]
    out$var <- NULL
    
  }
  
  return(out)
}
