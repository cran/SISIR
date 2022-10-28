#' @import magrittr
#' @importFrom dplyr mutate filter arrange select summarise ungroup inner_join
#' @importFrom tidyselect starts_with
#' @importFrom tidyr gather spread
#' @importFrom purrr map_dfc
#' @importFrom mixOmics pls
#' @importFrom stats sd as.hclust
#' @importFrom rlang .data

summary_pls <- function(dataset, target, groups) {
  # Description #####
  # returns a list of datasets with the same rows than 'dataset' and with a
  # number of columns corresponding to the different values of groups
  # (containing the PLS with respect to the target of the columns in 'dataset'
  # that are in the same group)
  # Inputs #####
  # dataset: a data.frame with observations in rows and variables in columns
  # target: target variable to predict
  # groups: a list with a hierarchy of groups
  # Outputs #####
  # a list where each entry is a dataset where the summary has been performed
  # for the groups
  
  p <- ncol(dataset)
  colnames(dataset) <- paste0("Var", 1:p)
  
  summary_data <- lapply(groups, function(cur_clust) {
    cur_data <- dataset %>% dplyr::mutate("sim" = 1:nrow(dataset)) %>% 
      tidyr::gather("Var", "value", - .data$sim)
    var_number <- as.numeric(gsub("Var", "", cur_data$Var))
    cur_data$clust <- cur_clust[var_number]
    
    pls_summary <- purrr::map_dfc(unique(cur_clust), compute_pls_summary, 
                                  dataset = cur_data, target = target)
    
    return(pls_summary)
  })
  
  return(summary_data)
}

compute_pls_summary <- function(aclust, dataset, target) {
  tmp <- dataset %>% dplyr::filter(.data$clust == aclust) %>% 
    tidyr::spread(.data$Var, .data$value) %>% dplyr::arrange(.data$sim) %>% 
    dplyr::select(tidyselect::starts_with("Var"))
  if (length(unique(dataset$Var[dataset$clust == aclust])) > 1) {
    out <- mixOmics::pls(tmp, target, ncomp = 1)
    return(unlist(out$variates$X))
  } 
  
  return(unlist(tmp))
}

summary_basics <- function(dataset, groups) {
  # Description #####
  # returns a list of datasets with the same rows than 'dataset' and with a
  # number of columns corresponding to the different values of groups
  # (containing the mean and sd of the columns within the same group)
  # Inputs #####
  # dataset: a data.frame with observations in rows and variables in columns
  # groups: a list with a hierarchy of groups
  # Outputs #####
  # a list where each entry is a dataset where the summary has been performed
  # for the groups
  
  p <- ncol(dataset)
  colnames(dataset) <- paste0("Var", 1:p)
  
  summary_data <- lapply(groups, function(cur_clust) {
    cur_data <- dataset %>% dplyr::mutate("sim" = 1:nrow(dataset)) %>% 
      tidyr::gather("Var", "value", - .data$sim)
    var_number <- as.numeric(gsub("Var", "", cur_data$Var))
    cur_data$clust <- cur_clust[var_number]
    
    mean_summary <- cur_data %>% dplyr::group_by(.data$sim, .data$clust) %>%
      dplyr::summarise(VarMean = mean(.data$value)) %>%
      tidyr::spread(.data$clust, value = .data$VarMean, sep = "mean") %>%
      dplyr::ungroup()
    
    sd_summary <- cur_data %>% dplyr::group_by(.data$sim, .data$clust) %>%
      dplyr::summarise(VarSD = sd(.data$value)) %>%
      tidyr::spread(.data$clust, value = .data$VarSD, sep = "sd") %>%
      dplyr::ungroup()
    
    cur_data <- dplyr::inner_join(mean_summary, sd_summary, by = "sim") %>%
      dplyr::select(- .data$sim)
    
    # remove column with NAs (for standard deviations)
    contains_na <- colSums(apply(cur_data, 2, is.na))
    cur_data <- cur_data[, contains_na == 0]
    
    return(cur_data)
  })
  
  return(summary_data)
}

summary_cclustofvar <- function(dataset, groups) {
  # Description #####
  # returns a list of datasets with the same rows than 'dataset' and with a
  # number of columns corresponding to the different values of groups
  # (containing the cov consensus variable)
  # Inputs #####
  # dataset: a data.frame with observations in rows and variables in columns
  # groups: a list with a hierarchy of groups
  # Outputs #####
  # a list where each entry is a dataset where the summary has been performed
  # for the groups
  
  p <- ncol(dataset)
  colnames(dataset) <- paste0("Var", 1:p)
  
  summary_data <- lapply(groups, function(cur_clust) {
    cur_data <- dataset %>% dplyr::mutate("sim" = 1:nrow(dataset)) %>% 
      tidyr::gather("Var", "value", - .data$sim)
    var_number <- as.numeric(gsub("Var", "", cur_data$Var))
    cur_data$clust <- cur_clust[var_number]
    
    cov_summary <- purrr::map_dfc(unique(cur_clust), compute_cov_summary, 
                                  dataset = cur_data)
    
    return(cov_summary)
  })
  
  return(summary_data)
}

compute_cov_summary <- function(aclust, dataset) {
  tmp <- dataset %>% dplyr::filter(.data$clust == aclust) %>% 
    tidyr::spread(.data$Var, .data$value) %>% dplyr::arrange(.data$sim) %>% 
    dplyr::select(starts_with("Var"))
  if (length(unique(dataset$Var[dataset$clust == aclust])) > 1) {
    out <- svd(tmp, 1, 1)
    n <- nrow(tmp)
    out <- sqrt(n) * out$d[1] * out$u[ ,1]
    return(out)
  }
  
  return(unlist(tmp))
}
