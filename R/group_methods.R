#' @importFrom stats cor
#' @importFrom adjclust adjClust cutree_chac

group_adjclust <- function(dataset) {
  # Description #####
  # returns a hierarchy of groups (list of clusterings) obtained by cutting the 
  # dendrogram provided by adjclust on 'dataset'
  # Inputs #####
  # dataset: a data.frame with observations in rows and variables in columns
  # Outputs #####
  # a list where each entry is a clustering named by the number of clusters
  # ordered from the largest number of clusters (number of variables) to 1
  
  corr_inputs <- stats::cor(dataset)
  hclust_res <- adjclust::adjClust(corr_inputs, type = "similarity")
  groups <- lapply(1:ncol(dataset), function(k) 
    adjclust::cutree_chac(hclust_res, k = k))
  hclust_res <- as.hclust(hclust_res)
  names(groups) <- 1:ncol(dataset)
  groups <- groups[ncol(dataset):1]
  
  return(list("groups" = groups, "dendro" = hclust_res))
}

group_cclustofvar <- function(dataset) {
  # Description #####
  # returns a hierarchy of groups (list of clusterings) obtained by cutting the 
  # dendrogram provided by an order constrained version of clustofvar on 
  # 'dataset'
  # Inputs #####
  # dataset: a data.frame with observations in rows and variables in columns
  # Outputs #####
  # a list where each entry is a clustering named by the number of clusters
  # ordered from the largest number of clusters (number of variables) to 1
  
  # initialization
  p <- ncol(dataset)
  inputs <- scale(dataset) / sqrt(nrow(dataset)-1)
  clust_res <- data.frame("number" = rep(0, p-1),
                          "left" = rep(0, p-1),
                          "right" = rep(0, p-1),
                          "link" = rep(0, p-1))
  
  all_groups <- as.list(1:p)
  names(all_groups) <- as.character(1:p)
  old_group_homogeneity <- rep(1, length(all_groups)) # because initial homogeneity is 1
  groups_homogeneity <- sapply(1:(length(all_groups)-1), function(ind) {
    selected <- c(all_groups[[ind]], all_groups[[ind+1]])
    homogeneity(dataset[, selected])
  })
  
  # loop over the number of clusters for merges
  for (level in 1:(p-1)) {
    ng <- length(old_group_homogeneity)
    criteria <- old_group_homogeneity[1:(ng-1)] + old_group_homogeneity[2:ng] - 
      groups_homogeneity 
    to_group <- which.min(criteria)
    cur_nei <- as.numeric(names(all_groups)[c(to_group, to_group+1)])
    clust_res[level, ] <- c(-level, cur_nei, min(criteria))
    # update left and right homogeneity
    old_group_homogeneity <- groups_homogeneity
    if (level < p-1) {
      if (to_group > 1) {
        selected <- c(all_groups[[(to_group-1)]], all_groups[[to_group]],
                      all_groups[[to_group+1]])
        left_homo <- homogeneity(dataset[, selected])
        groups_homogeneity[to_group - 1] <- left_homo
        if (to_group < length(all_groups) - 1) {
          selected <- c(all_groups[[(to_group+2)]], all_groups[[to_group]],
                        all_groups[[to_group+1]])
          right_homo <- homogeneity(dataset[, selected])
          groups_homogeneity[to_group] <- right_homo
          groups_homogeneity <- groups_homogeneity[-c(to_group+1)]
        } else {
          groups_homogeneity <- groups_homogeneity[-to_group]
        }
      } else {
        selected <- c(all_groups[[(to_group+2)]], all_groups[[to_group]],
                      all_groups[[to_group+1]])
        right_homo <- homogeneity(dataset[, selected])
        groups_homogeneity[to_group+1] <- right_homo
        groups_homogeneity <- groups_homogeneity[-1]
      }
    }
    
    # update groups
    all_groups[[to_group]] <- c(all_groups[[to_group]], 
                                all_groups[[to_group+1]])
    names(all_groups)[to_group] <- -level
    all_groups <- all_groups[-c(to_group+1)]
  }
  
  clust_res[, 2:3] <- - clust_res[, 2:3]
  out_clust <- list("merge" = matrix(c(as.integer(clust_res[, 2]),
                                       as.integer(clust_res[, 3])),
                                     ncol = 2, byrow = FALSE),
                    "height" = clust_res$link,
                    "order" = 1:p,
                    "labels" = names(dataset),
                    "method" = "CClustOfVar",
                    "call" = NULL,
                    "dist.method" = "PCAmix")
  class(out_clust) <- "chac"
  
  groups <- lapply(1:ncol(dataset), function(k) 
    adjclust::cutree_chac(out_clust, k = k))
  names(groups) <- 1:ncol(dataset)
  groups <- groups[ncol(dataset):1]
  
  out_clust <- as.hclust(out_clust)
  
  return(list("groups" = groups, "dendro" = out_clust))
}

# useful functions #####
homogeneity <- function(mat) { 
  # compute homogeneity criterion (first svd value)
  out <- svd(mat, 1, 1)
  return(out$d[1])
}