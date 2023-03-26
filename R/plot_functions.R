#' @import ggplot2
#' @importFrom reshape2 melt
#' @import dendextend
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats as.dendrogram na.omit
#' @importFrom graphics abline

plot_dendrogram <- function(x) {
  
  dendro <- as.dendrogram(x$dendro)
  nlev <- dendextend::nleaves(dendro)
  labs <- labels(dendro)
  dendro <- set_labels(dendro, rep("", length(labs)))
  dendro %>% plot(leaflab = "none") 
  
  all_at <- unique(x$mse$clust)
  nrange <- length(x$groups)
  kmax <- max(all_at)
  basics <- ncol(x$importances[[1]]) != (nrow(x$mse) / length(all_at))
  if (basics) {
    all_groups <- rep(names(x$groups), 2)
    all_types <- rep(c("mean", "sd"), each = nrange)
  }
  
  all_h <- rev(x$dendro$height)
  
  # rectangles or lines
  if (all(diff(x$dendro$height) >= 0)) {
    dendro %>% rect.dendrogram(k = kmax)
  } else {
    message(paste("Reversals detected in the dendrogram. Rectangles are not",
                  "relevant and thus, they are not displayed."))
  }
  
  # horizontal lines ('range_at')
  if (nrange > 1) {
    for (k in setdiff(all_at, max(all_at))) {
      abline(h = (all_h[k] + all_h[k-1]) / 2, col = "darkred", lwd = 2)
    }
  }
  
  if (!basics) {
    if (nrange > 3) {
      cutb <- 3
      message(paste("Only the first 3 selections are represented below. If ",
                    "you want to display others, please use `extract_at`",
                    "before plotting."))
    } else cutb <- nrange
  } else {
    if (nrange > 2) {
      cutb <- 2
      message(paste("Only the first 2 selections are represented below. If",
                    "you want to display others, please use `extract_at`",
                    "before plotting."))
    } else cutb <- nrange
    all_types_c <- rep(c("mean", "sd"), each = cutb)
    all_groups_c <- rep(names(x$groups)[1:cutb], 2)
  }
  
  # bars
  if ("selected" %in% names(x)) { # case 'selected'
    values <- replicate(cutb, rep("grey", length(labs)), simplify = FALSE)
    if (!basics) {
      values <- mapply(function(a, b) {
        names(a) <- labs; a[b] <- "darkred"; return(a)
        }, values, x$selected[1:cutb])
      bar_labs <- paste("#int.:", names(x$groups)[1:cutb], "- sel:")
    } else {
      values_mean <- mapply(function(a, b) {
        names(a) <- labs; a[b$variable[b$"in.mean"]] <- "darkred"; return(a)
      }, values, x$selected[1:cutb])
      values_sd <- mapply(function(a, b) {
        names(a) <- labs; a[b$variable[b$"in.sd"]] <- "darkred"; return(a)
      }, values, x$selected[1:cutb])
      values <- cbind(values_mean, values_sd)
      bar_labs <- paste(all_groups_c, "/", all_types_c, "sel.")
    }
  } else {
    best <- tapply(x$mse$mse, x$mse$clust, which.min)
    if (!basics) {
      values <- mapply(function(a, b) return(a[, b]), 
                       x$importances[1:cutb], best[1:cutb])
      values <- apply(values, 2, cut, breaks = 9, labels = FALSE)
      values <- apply(values, 2, function(acol) brewer.pal(9, "YlOrRd")[acol])
      bar_labs <- paste("#int.:", names(x$groups)[1:cutb], "- imp.")
    } else {
      bestm <- paste0("mean", best)
      bestsd <- paste0("sd", best)
      values_m <- mapply(function(a, b) return(a[, b]), 
                         x$importances[1:cutb], bestm[1:cutb])
      values_m <- apply(values_m, 2, cut, breaks = 9, labels = FALSE)
      values_m <- apply(values_m, 2, 
                        function(acol) brewer.pal(9, "YlOrRd")[acol])
      values_sd <- mapply(function(a, b) return(a[, b]), 
                          x$importances[1:cutb], bestsd[1:cutb])
      values_sd <- apply(values_sd, 2, cut, breaks = 9, labels = FALSE)
      values_sd <- apply(values_sd, 2, 
                         function(acol) brewer.pal(9, "YlOrRd")[acol])
      values_sd[is.na(values_sd)] <- "black"
      values <- cbind(values_m, values_sd)
      bar_labs <- paste(all_groups_c, "/", all_types_c, "- imp.")
    }
  }
  # bars
  colored_bars(colors = values, dend = dendro, rowLabels = bar_labs)

  return(invisible(NULL))
}

plot_importance <- function(x, shape.imp) {
  lobj <- length(x$importances)
  nrepeats <- nrow(x$mse) / lobj
  basics <- ncol(x$importances[[1]]) != nrepeats
  
  suppressMessages({ df <- reshape2::melt(x$importances) })
  if (!basics) { # not basics
    names(df) <- c("variable", "repeats", "importance", "at")
  } else { # basics
    df$type <- rep("mean", nrow(df))
    df$type[grep("sd", df$variable)] <- "sd"
    names(df) <- c("repeats", "importance", "at", "type")
  }
  
  if (shape.imp == "boxplot") {
    p <- ggplot(df, aes(x = .data$repeats, y = .data$importance)) + 
      geom_boxplot() + theme_bw()
  } else if (shape.imp == "histogram") {
    p <- ggplot(df, aes(x = .data$importance)) + geom_histogram(bins = 30) + 
      theme_bw()
  }
  
  if (lobj > 1) {
    if (!basics | shape.imp == "boxplot") {
      p <- p + facet_wrap(~ at)
    } else {
      p <- p + facet_grid(type ~ at)
    }
  } else if (shape.imp == "histogram" & basics) { # basics
    p <- p + facet_wrap(~ type)
  }
  
  return(p)
}

plot_selection <- function(x, sel.type, threshold) {
  is_selected <- "selected" %in% names(x)
  
  if (threshold != "none") {
    if (!is.numeric(threshold)) 
      stop("'threshold' must be numeric.", call. = FALSE)
    
    if (is_selected) 
      stop(paste("A selection method has already been used: 'threshold' must",
                 "be 'none'."),
           call. = FALSE)
           
    if (sel.type != "selection") {
      warning(paste("A 'threshold' has passed to the plot function while",
                    "'sel.type' is not 'selection'. Automatically switching it",
                    "to 'selection'."),
              call. = FALSE)
      sel.type <- "selection"
    }
  }
  
  if (sel.type == "selection" && !is_selected && threshold == "none")
    stop(paste("Variable selection was not used on this result. Choose",
               "sel.type = 'importance'"),
         call. = FALSE)
  
  lobj <- length(x$importances)
  nrepeats <- nrow(x$mse) / lobj
  basics <- ncol(x$importances[[1]]) != nrepeats
  var_names <- names(x$groups[[1]])
  
  if (threshold != "none") {
    if (!basics) { # not basics, 'at'
      best <- tapply(x$mse$mse, x$mse$clust, which.min)
      x$selected <- mapply(function(a, b) { rownames(a)[a[, b] > threshold] }, 
                           x$importances, best, SIMPLIFY = FALSE)
    } else {
      best <- tapply(x$mse$mse, x$mse$clust, which.min)
      x$selected <- mapply(function(a, b) {
        bmean <- paste0("mean", b)
        bmean <- rownames(a)[a[, bmean] > threshold]
        bsd <- paste0("sd", b)
        bsd <- rownames(a)[a[, bsd] > threshold]
        bsd <- na.omit(bsd)
        variables <- unique(c(bmean, bsd))
        variables <- var_names[var_names %in% variables]
        out <- data.frame("variable" = variables,
                          "in.mean" = variables %in% bmean,
                          "in.sd" = variables %in% bsd)
        return(out)
      }, x$importances, best, SIMPLIFY = FALSE)
    }
  }
  
  if (sel.type == "selection") { # not basics, 'at'
    if (lobj == 1 & !basics) {
      dfline <- make_line_df(x$selected[[1]], var_names)
      dfrect <- make_rect_df(x$selected[[1]], var_names)
      
      p <- ggplot() +
        geom_rect(data = dfrect, 
                  aes(xmin = .data$xstart, ymin = .data$ystart, 
                      xmax = .data$xend, ymax = .data$yend), 
                  alpha = 0.2)
    } else if (!basics) { # not basics, 'range_at'
      all_levels <- seq(0.7, 1, length.out = length(x$selected))
      all_levels2 <- seq(0, 0.3, length.out = length(x$selected))
      dfline <- mapply(function(a, b, d) {
        make_line_df(a, var_names = var_names, level = b, level2 = d)
      }, x$selected, all_levels, all_levels2, SIMPLIFY = FALSE)
      
      all_at <- sapply(dfline, nrow)
      all_at <- rep(names(dfline), all_at)
      dfline <- data.frame(Reduce(rbind, dfline), "at" = all_at)
      
      p <- ggplot(data = dfline, aes(group = .data$at, colour = .data$at))
    } else if (lobj == 1) { # basics, 'at'
      sel_mean <- x$selected[[1]]$variable[x$selected[[1]]$"in.mean"]
      dfline_mean <- make_line_df(sel_mean, var_names, level2 = 0.2)
      sel_sd <- x$selected[[1]]$variable[x$selected[[1]]$"in.sd"]
      dfline_sd <- make_line_df(sel_sd, var_names, level = 0.8)
      dfline <- rbind(dfline_mean, dfline_sd)
      dfline$type <- c(rep("mean", nrow(dfline_mean)),
                       rep("sd", nrow(dfline_sd)))
      
      p <- ggplot(data = dfline, aes(group = .data$type, colour = .data$type))
    } else { # basics, 'range_at'
      dfmean <- lapply(x$selected, 
                       function(alist) alist$variable[alist$"in.mean"])
      names(dfmean) <- paste0("mean", names(dfmean))
      dfsd <- lapply(x$selected, function(alist) alist$variable[alist$"in.sd"])
      names(dfsd) <- paste0("sd", names(dfmean))
      df <- c(dfmean, dfsd)
      
      all_levels <- seq(0.7, 1, length.out = length(df))
      all_levels2 <- seq(0, 0.3, length.out = length(df))
      dfline <- mapply(function(a, b, d) {
        make_line_df(a, var_names = var_names, level = b, level2 = d)
      }, df, all_levels, all_levels2, SIMPLIFY = FALSE)
      
      all_at <- sapply(dfline, nrow)
      all_at <- rep(names(dfline), all_at)
      dfline <- data.frame(Reduce(rbind, dfline), "at" = all_at)
      
      p <- ggplot(data = dfline, aes(group = .data$at, colour = .data$at))
    }
    
     p <- p +
      geom_segment(data = dfline, 
                   aes(x = .data$xstart, y = .data$ystart, xend = .data$xend, 
                       yend = .data$yend), 
                   size = 1) +
      theme_bw() + xlim(1, length(var_names)) + 
      scale_y_continuous(breaks = c(0, 1), minor_breaks = c(0, 1), 
                         labels = c("not selected", "selected"),
                         limits = c(0, 1.23)) +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank())
     if (lobj == 1 & basics) { # basics, 'at'
       p <- p + guides(colour = guide_legend(title = "for summary..."))
     } else if (!basics & lobj > 1) { # not basics, 'range_at'
       p <- p + guides(colour = guide_legend(title = "# intervals"))
     } else if (lobj > 1) { # basics, 'range_at'
       p <- p + guides(colour = guide_legend(title = "summary + # interv."))
     }
    
  } else { # sel.type == "importance"
    if (lobj == 1 & !basics) {
      imp_names <- rownames(x$importances[[1]])
      impmean <- data.frame("x" = match(imp_names, var_names),
                            "mean" = apply(x$importances[[1]], 1, mean))
      impmin <- data.frame("x" = match(imp_names, var_names),
                           "min" = apply(x$importances[[1]], 1, min))
      impmax <- data.frame("x" = match(imp_names, var_names),
                           "max" = apply(x$importances[[1]], 1, max))
      df <- merge(impmean, impmin)
      df <- merge(df, impmax)
      df <- fill_na(df, c("mean", "min", "max"))
      p <- ggplot(df, aes(x = .data$x))
    } else if (!basics) {
      impmean <- compute_sumimp(x$importances, "mean", var_names)
      impmin <- compute_sumimp(x$importances, "min", var_names)
      impmax <- compute_sumimp(x$importances, "max", var_names)
      df <- merge(impmean, impmin)
      df <- merge(df, impmax)
      df <- fill_na(df, c("mean", "min", "max"))
      all_at <- as.character(sort(unique(df$at)))
      df$at <- factor(df$at, levels = all_at, ordered = TRUE)
      p <- ggplot(df, aes(x = .data$x, group = .data$at, colour = .data$at, 
                          fill = .data$at))
    } else if (lobj == 1) {
      impmean <- compute_sumimp1_basics(x$importances[[1]], "mean", var_names)
      impmin <- compute_sumimp1_basics(x$importances[[1]], "min", var_names)
      impmax <- compute_sumimp1_basics(x$importances[[1]], "max", var_names)
      df <- merge(impmean, impmin)
      df <- merge(df, impmax)
      df <- fill_na(df, c("mean", "min", "max"))
      p <- ggplot(df, aes(x = .data$x, colour = .data$type, fill = .data$type, 
                          group = .data$type))
    } else {
      impmean <- compute_sumimp_basics(x$importances, "mean",var_names)
      impmin <- compute_sumimp_basics(x$importances, "min", var_names)
      impmax <- compute_sumimp_basics(x$importances, "max", var_names)
      df <- merge(impmean, impmin)
      df <- merge(df, impmax)
      df <- fill_na(df, c("mean", "min", "max"))
      p <- ggplot(df, aes(x = .data$x, group = .data$at, colour = .data$at, 
                          fill = .data$at))
    }
    
    p <- p + geom_ribbon(aes(ymin = .data$min, ymax = .data$max), alpha = 0.1) + 
      geom_line(aes(y = mean)) + theme_bw() + 
      ylab("variable importance") +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  return(p)
}

plot_quality <- function(x, quality.crit) {
  
  valid_criteria <- "mse"
  if ("computational.times" %in% names(x))  {
    valid_criteria <- c(valid_criteria, "time")
  }
  if ("quality" %in% names(x)) {
    valid_criteria <- c(valid_criteria, "Precision", "Recall", "ARI", "NMI")
  }
  
  crit_ok <- all(sapply(quality.crit, function(cc) cc %in% valid_criteria))
  if (!crit_ok || length(quality.crit) > 2) {
    stop(paste0("'quality.crit' must be a vector with length at most 2 in ",
                paste(valid_criteria, collapse = ", "), "."),
         call. = FALSE)
  }
  
  if (length(quality.crit) == 2 && "time" %in% quality.crit) {
    stop("'time' is a valid quality criterion to plot only taken alone.",
         call. = FALSE)
  }
  
  if ("time" %in% quality.crit) {
    df <- data.frame("criterion" = c(unname(x$"computational.times")),
                     "step" = names(x$"computational.times"))
    df$step <- factor(df$step, levels = names(x$"computational.times"),
                      ordered = TRUE)
    p <- ggplot(df, aes(fill = .data$step, y = .data$criterion, x = 1)) + 
      geom_bar(stat = "identity") + theme_bw() + xlim(0, 2) +
      ylab("computational time (s)") +
      theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) +
      scale_fill_brewer(type = "qual", palette = 7)
    return(p)
  }
  
  if (length(quality.crit) == 1) {
    # build data frame
    if (quality.crit == "mse") {
      df <- data.frame("criterion" = x$mse$mse, "at" = x$mse$clust)
      ylimits <- c(0, max(df$criterion))
    } else {
      df <- data.frame("criterion" = x$quality[, quality.crit], 
                       "at" = as.numeric(x$quality$clust))
      if (quality.crit == "ARI") {
        ylimits <- c(-1, 1)
      } else ylimits <- c(0, 1)
    }
    
    # make plot
    p <- ggplot(df, aes(x = .data$at, y = .data$criterion)) + 
      geom_jitter(width=0.2, height = 0) + theme_bw() + 
      xlab("number of intervals") + ylab(quality.crit) + ylim(ylimits) +
      scale_x_continuous(breaks = unique(df$at), 
                         limits = c(min(df$at) - 0.5, max(df$at + 0.5)))
    
  } else {
    if ("mse" %in% quality.crit) {
      quality.crit <- setdiff(quality.crit, "mse")
      df <- data.frame("x" = x$quality[, quality.crit], "y" = x$mse$mse, 
                       "at" = as.factor(x$quality$clust))
      quality_names <- c(quality.crit, "mse")
      if (quality.crit == "ARI") {
        xlimits <- c(-1, 1)
      } else xlimits <- c(0, 1)
      ylimits <- c(0, max(df$y))
    } else {
      df <- data.frame("x" = x$quality[, quality.crit[1]], 
                       "y" = x$quality[, quality.crit[2]], 
                       "at" = as.factor(x$quality$clust))
      quality_names <- quality.crit
      if (quality.crit[1] == "ARI") {
        xlimits <- c(-1, 1)
      } else xlimits <- c(0, 1)
      if (quality.crit[2] == "ARI") {
        ylimits <- c(-1, 1)
      } else ylimits <- c(0, 1)
    }
    p <- ggplot(df, aes(x = .data$x, y = .data$y, colour = .data$at)) + 
      geom_point() + theme_bw() + xlab(quality_names[1]) + 
      ylab(quality_names[2]) + xlim(xlimits) + ylim(ylimits) + 
      scale_colour_discrete(name = "# intervals")
    if (length(unique(df$at)) > 20) {
      p <- p + theme(legend.position = "none")
    }
  }
  
  return(p)
}

compute_sumimp <- function(importances, func_name, var_names) {
  FUN <- eval(func_name)
  impsummary <- lapply(importances, function(alist) apply(alist, 1, FUN))
  impsummary <- lapply(impsummary, function(alist) {
    out <- rep(NA, length(var_names))
    names(out) <- var_names
    out[match(names(alist), var_names)] <- alist
    return(out)
  })
  impsummary <- Reduce(cbind, impsummary)
  rownames(impsummary) <- match(rownames(impsummary), var_names)
  colnames(impsummary) <- names(importances)
  impsummary <- melt(impsummary)
  names(impsummary) <- c("x", "at", func_name)
  return(impsummary)
}

compute_sumimp1_basics <- function(importances, func_name, var_names) {
  FUN <- eval(func_name)
  
  selected <- grep("mean", names(importances))
  importances_mean <- importances[, selected]
  selected <- grep("sd", names(importances))
  importances_sd <- importances[, selected]
  importances <- list(importances_mean, importances_sd)
  names(importances) <- c("mean", "sd")
  
  impsummary <- lapply(importances, function(alist) 
    data.frame("x" = match(rownames(alist), var_names),
               "importance" = apply(alist, 1, FUN, na.rm = TRUE)))
  impsummary <- reshape2::melt(impsummary, id = "x")
  names(impsummary) <- c("x", "variable", func_name, "type")
  
  return(impsummary)
}

compute_sumimp_basics <- function(importances, func_name, var_names) {
  FUN <- eval(func_name)
  importances_mean <- lapply(importances, function(alist) {
    selected <- grep("mean", names(alist))
    return(alist[, selected])
  })
  importances_sd <- lapply(importances, function(alist) {
    selected <- grep("sd", names(alist))
    return(alist[, selected])
  })
  importances <- c(importances_mean, importances_sd)
  imp_names <- rep(c("mean", "sd"), each = length(importances) / 2)
  names(importances) <- paste0(imp_names, names(importances))
  impsummary <- lapply(importances, function(alist) 
    apply(alist, 1, FUN, na.rm = TRUE))
  impsummary <- lapply(impsummary, function(alist) {
    out <- rep(NA, length(var_names))
    names(out) <- var_names
    out[match(names(alist), var_names)] <- alist
    return(out)
  })
  impsummary <- Reduce(cbind, impsummary)
  rownames(impsummary) <- match(rownames(impsummary), var_names)
  colnames(impsummary) <- names(importances)
  impsummary <- melt(impsummary)
  names(impsummary) <- c("x", "at", func_name)
  return(impsummary)
}

fill_na <- function(dataset, var) {
  out <- sapply(seq_along(var), function(ind) {
    outcol <- dataset[, var[ind]]
    outcol[is.na(outcol)] <- 0
    return(outcol)
  })
  dataset[, var] <- out
  out <- sapply(seq_along(var), function(ind) {
    outcol <- dataset[, var[ind]]
    outcol[outcol == Inf] <- 0
    return(outcol)
  })
  dataset[, var] <- out
  return(dataset)
}

make_line_df <- function(selected, var_names, level = 1, level2 = 0) {
  which_sel <- match(selected, var_names)
  if (length(which_sel) > 0 & length(which_sel) < length(var_names)) {
    diff_which <- diff(which_sel)
    starts <- which(diff_which > 1) + 1
    starts <- c(which_sel[1], which_sel[starts])
    ends <- which(diff_which > 1)
    ends <- c(which_sel[ends], which_sel[length(which_sel)])
    df <- data.frame("xstart" = starts, "ystart" = rep(level, length(starts)),
                     "xend" = ends, "yend" = rep(level, length(starts)))
    which_sel <- setdiff(seq_along(var_names), which_sel)
    diff_which <- diff(which_sel)
    starts <- which(diff_which > 1) + 1
    starts <- c(which_sel[1], which_sel[starts])
    ends <- which(diff_which > 1)
    ends <- c(which_sel[ends], which_sel[length(which_sel)])
    df2 <- data.frame("xstart" = starts, "ystart" = rep(level2, length(starts)),
                      "xend" = ends, "yend" = rep(level2, length(starts)))
    df <- rbind(df, df2)
  } else if (length(which_sel) == 0) {
    df <- data.frame("xstart" = 1, "ystart" = 0, "xend" = length(var_names),
                     "yend" = 0)
  } else {
    df <- data.frame("xstart" = 1, "ystart" = level, "xend" = length(var_names),
                     "yend" = level)
  }
  return(df)
}

make_rect_df <- function(selected, var_names) {
  which_sel <- match(selected, var_names)
  diff_which <- diff(which_sel)
  starts <- which(diff_which > 1) + 1
  starts <- c(which_sel[1], which_sel[starts])
  ends <- which(diff_which > 1)
  ends <- c(which_sel[ends], which_sel[length(which_sel)])
  df_rect <- data.frame("xstart" = starts, "ystart" = rep(0, length(starts)), 
                        "xend" = ends, "yend" = rep(1, length(starts)))

  return(df_rect)
}
