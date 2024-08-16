################################################################################
# sparse SIR
################################################################################
#' @import glmnet
#' @importFrom stats predict
#' 
#' @title sparse SIR
#' @export
#'
#' @description
#' \code{sparseSIR} performs the second step of the method (shrinkage of ridge
#' SIR results
#'
#' @param object an object of class \code{ridgeRes} as obtained from the 
#' function \code{\link{ridgeSIR}}
#' @param inter_len (numeric) vector with interval lengths
#' @param adaptive should the function returns the list of strong zeros and non
#' strong zeros (logical). Default to FALSE
#' @param sel_prop used only when \code{adaptive = TRUE}. Fraction of the 
#' coefficients that will be considered as strong zeros and strong non zeros.
#' Default to 0.05
#' @param parallel whether the computation should be performed in parallel or
#' not. Logical. Default is FALSE
#' @param ncores number of cores to use if \code{parallel = TRUE}. If left to 
#' NULL, all available cores minus one are used
#' 
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' 
#' @references {Picheny, V., Servien, R., and Villa-Vialaneix, N. (2019) 
#' Interpretable sparse SIR for digitized functional data. 
#' \emph{Statistics and Computing}, \strong{29}(2), 255--267.}
#' 
#' @seealso \code{\link{ridgeSIR}}, \code{\link{project.sparseRes}}, 
#' \code{\link{SISIR}}
#' 
#' @examples
#' set.seed(1140)
#' tsteps <- seq(0, 1, length = 200)
#' nsim <- 100
#' simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
#' x <- t(replicate(nsim, simulate_bm()))
#' beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
#' beta[((tsteps < 0.2) | (tsteps > 0.5)), 1] <- 0
#' beta[((tsteps < 0.6) | (tsteps > 0.75)), 2] <- 0
#' y <- log(abs(x %*% beta[ ,1]) + 1) + sqrt(abs(x %*% beta[ ,2]))
#' y <- y + rnorm(nsim, sd = 0.1)
#' res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = 10^8)
#' res_sparse <- sparseSIR(res_ridge, rep(10, 20))
#' 
#' @return S3 object of class \code{sparseRes}: a list consisting of
#' \describe{
#'    \item{\code{sEDR}}{ the estimated EDR space (a p x d matrix)}
#'    \item{\code{alpha}}{ the estimated shrinkage coefficients (a vector having
#'    a length similar to \code{inter_len})}
#'    \item{\code{quality}}{ a vector with various qualities for the model (see
#'    Details)}
#'    \item{\code{adapt_res}}{ if \code{adaptive = TRUE}, a list of two vectors: 
#'    \describe{
#'      \item{\code{nonzeros}}{ indexes of variables that are strong non zeros}
#'      \item{\code{zeros}}{ indexes of variables that are strong zeros}
#'    }}
#'    \item{\code{parameters}}{ a list of hyper-parameters for the method: 
#'    \describe{
#'      \item{\code{inter_len}}{ lengths of intervals}
#'      \item{\code{sel_prop}}{ if \code{adaptive = TRUE}, fraction of the 
#'      coefficients which are considered as strong zeros or strong non zeros}
#'    }}
#'    \item{\code{rSIR}}{ same as the input \code{object}}
#'    \item{\code{fit}}{ a list for LASSO fit with:
#'    \describe{
#'      \item{\code{glmnet}}{ result of the \code{\link[glmnet]{glmnet}} 
#'      function}
#'      \item{\code{lambda}}{ value of the best Lasso parameter by CV}
#'      \item{\code{x}}{ exploratory variable values as passed to fit the 
#'      model}
#'    }}
#'  }
#'  
#'  @details Different quality criteria used to select the best models among a
#'  list of models with different interval definitions. Quality criteria are:
#'  log-likelihood (\code{loglik}), cross-validation error as provided by the
#'  function \code{\link[glmnet]{glmnet}}, two versions of the AIC (\code{AIC} 
#'  and \code{AIC2}) and of the BIC (\code{BIC} and \code{BIC2}) in which the 
#'  number of parameters is either the number of non null intervals or the 
#'  number of non null parameters with respect to the original variables.

sparseSIR <- function(object, inter_len, adaptive = FALSE, sel_prop = 0.05,
                      parallel = FALSE, ncores = NULL) {
  if (parallel) {
    if (is.null(ncores)) ncores <- min(detectCores() - 1)
    print(ncores)
    registerDoParallel(cores = ncores)
  }
  oldwarn <- getOption("warn")
  options(warn = -1)
  
  H <- object$parameters$H
  D <- length(inter_len)
  p <- nrow(object$EDR)
  if (sum(inter_len) != p)
    stop("Sum of 'inter_len' must be equal to the number of columns in original dataset.")
  d <- ncol(object$EDR)
  cmean_x <- object$utils$cmean_x
  slices <- object$utils$slices
  rEDR <- object$EDR
  x <- object$data$x
  
  coef_elastic <- 1 # note: can be passed as an argument for elastinet
  
  # definition of the dependant variable: vector of length dn
  Proj <- crossprod(sweep(cmean_x[ ,slices], 1, apply(x, 2, mean), "-") , rEDR)
  Proj <- as.vector(Proj)
  
  # definition of the explanatory variables: matrix of dimension (dn)xD
  explanatory <- apply(rEDR, 2, function(acol) {
    Delta <- matrix(0, p, D)
    Delta[cbind(1:p, rep(1:D, inter_len))] <- acol
    return(list("out" = x %*% Delta))
  })
  explanatory <- lapply(explanatory, function(alist) alist$out)
  explanatory <- Reduce(rbind, explanatory)
  
  # model fit (CV)
  fit <- cv.glmnet(x = explanatory, y = Proj, alpha = coef_elastic, 
                   parallel = parallel)
  if (parallel) {
    stopImplicitCluster()
  }
  if (length(fit$lambda) < 20) {# in this case, fit on a very fine grid
    lambda_seq <- fit$lambda[1] * 10^(-seq(1, 10, length = 1000))
    fit <- cv.glmnet(x = explanatory, y = Proj, alpha = coef_elastic, 
                     lambda = lambda_seq, parallel = parallel)
    if (parallel) {
      stopImplicitCluster()
    }
  }
  best_lambda <- fit$lambda.min
  ind_best <- which(best_lambda == fit$lambda.min)
  
  # re-fit with a finer grid to obtain the whole regularization path (no CV)
  seq_lambda <- c(as.vector(outer(seq(1:990)/100, 10**(1:-10), "*")),
                  best_lambda)
  seq_lambda <- sort(unique(seq_lambda), decreasing = TRUE)
  ind_best <- which(seq_lambda == best_lambda)
  
  best_fit <- glmnet(x = explanatory, y = Proj, alpha = coef_elastic,
                     lambda = seq_lambda, family="gaussian")
  alpha <- best_fit$beta[ ,ind_best]

  sEDR <- sweep(rEDR, 1, rep(alpha, inter_len), "*")
  
  # quality criteria TODO: see which ones to keep
  nobs <- fit$glmnet.fit$nobs
  loglik <- log(mean((predict(fit, explanatory, s=best_lambda) - Proj)^2))
  AIC <- nobs * loglik + 2 * sum(alpha != 0)
  BIC <- nobs * loglik + sum(alpha != 0) * log(nobs)
  nb_param <- sum(as.numeric(alpha != 0) * inter_len)
  AIC2 <- nobs * loglik + 2 * nb_param
  BIC2 <- nobs * loglik + nb_param * log(nobs)
  CVerror <- min(fit$cvm)
  quality <- c(loglik, AIC, BIC, AIC2, BIC2, CVerror)
  names(quality) <- c("loglik", "AIC", "BIC", "AIC2", "BIC2", "CVerror")
  
  if (adaptive) {# find hard/soft zeros
    ind_nonzero <- which.min(abs(best_fit$df - sel_prop*D))
    alpha_nonzero <- best_fit$beta[ ,ind_nonzero]
    ind_zero <- which.min(abs(best_fit$df - (1-sel_prop)*D))
    alpha_zero <- best_fit$beta[ ,ind_zero]
    adapt <- list("nonzeros" = alpha_nonzero, "zeros" = alpha_zero)
  } else adapt <- NULL
  
  parameters <- list("inter_len" = inter_len, "sel_prop" = sel_prop)

  res <- list("sEDR" = sEDR, "alpha" = unname(alpha), "quality" = quality, 
              "adapt_res" = adapt, "parameters" = parameters, "rSIR" = object,
              "fit" = list("glmnet" = best_fit, "lambda" = best_lambda,
                           "x" = explanatory))
  class(res) <- "sparseRes"
  
  options(warn = oldwarn)
  
  return(res)
}

################################################################################
# Methods for objects of class sparseRes
################################################################################
#' @title Print sparseRes object
#' @name sparseRes
#' @export
#' @aliases summary.sparseRes
#' @aliases summary.sparseRes
#' @aliases print.sparseRes
#' @aliases sparseRes-class
#' @description Print a summary of the result of \code{\link{sparseSIR}} (
#' \code{sparseRes} object)
#' @param object a \code{sparseRes} object
#' @param x a \code{sparseRes} object
#' @param ... not used
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inra.fr}}
#' @seealso \code{\link{sparseSIR}}
#' 
summary.sparseRes <- function(object, ...) {
  cat("Sparse SIR results with:\n\n",
      object$rSIR$parameters$H, "slices\n",
      "dimension of the EDR space is:", object$rSIR$parameters$d, "\n",
      "regularization parameter is:", object$rSIR$parameters$mu2, "\n\n")
  cat(sum(object$alpha != 0), "non zero coefficient(s) out of", 
      length(object$alpha), "possible.\n")
  cat("Shrinkage coefficients are in $alpha.", min(c(length(object$alpha), 10)),
      "first coefficients are:\n")
  print(head(object$alpha, 10))
  cat("\nThe EDR space is in '$sEDR'. First 10 rows are:\n\n")
  print(head(object$sEDR, 10))
}

#' @export
#' @rdname sparseRes
print.sparseRes <- function(x, ...) {
  summary.sparseRes(x)
}

#' @title sparse SIR
#' @name project
#' @export
#' @aliases project.sparseRes
#'
#' @description
#' \code{project} performs the projection on the sparse EDR space (as obtained
#' by the \code{\link[glmnet]{glmnet}})
#'
#' @param object an object of class \code{sparseRes} as obtained from the 
#' function \code{\link{sparseSIR}}
#' 
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' 
#' @references {Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016) 
#' Interpretable sparse SIR for digitized functional data.
#' \emph{Statistics and Computing}, \strong{29}(2), 255--267.}
#' 
#' @seealso \code{\link{sparseSIR}}
#' 
#' @details The projection is obtained by the function 
#' \code{\link[glmnet]{predict.glmnet}}.
#' 
#' @examples
#' set.seed(1140)
#' tsteps <- seq(0, 1, length = 200)
#' nsim <- 100
#' simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
#' x <- t(replicate(nsim, simulate_bm()))
#' beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
#' beta[((tsteps < 0.2) | (tsteps > 0.5)), 1] <- 0
#' beta[((tsteps < 0.6) | (tsteps > 0.75)), 2] <- 0
#' y <- log(abs(x %*% beta[ ,1]) + 1) + sqrt(abs(x %*% beta[ ,2]))
#' y <- y + rnorm(nsim, sd = 0.1)
#' \donttest{
#' res_ridge <- ridgeSIR(x, y, H = 10, d = 2)
#' res_sparse <- sparseSIR(res_ridge, rep(1, ncol(x)))
#' proj_data <- project(res_sparse)
#' }
#' 
#' @return a matrix of dimension n x d with the projection of the observations
#' on the d dimensions of the sparse EDR space
#' 
project.sparseRes <- function(object) {
  res <- predict.glmnet(object$fit$glmnet, object$fit$x, object$fit$lambda)
  res <- matrix(res, ncol = object$rSIR$parameters$d, byrow = FALSE)
  return(res)
}

#' @export
#' @rdname project
project <- function(object) {
  UseMethod("project")
}

################################################################################
# Iterative interval search
################################################################################
#' @title Interval Sparse SIR
#' @export
#'
#' @description
#' \code{SISIR} performs an automatic search of relevant intervals
#'
#' @param object an object of class \code{ridgeRes} as obtained from the 
#' function \code{\link{ridgeSIR}}
#' @param inter_len (numeric) vector with interval lengths for the initial 
#' state. Default is to set one interval for each variable (all intervals have
#' length 1)
#' @param sel_prop fraction of the coefficients that will be considered as 
#' strong zeros and strong non zeros. Default to 0.05
#' @param itermax maximum number of iterations. Default to Inf
#' @param minint minimum number of intervals. Default to 2
#' @param parallel whether the computation should be performed in parallel or
#' not. Logical. Default is FALSE
#' @param ncores number of cores to use if \code{parallel = TRUE}. If left to 
#' NULL, all available cores minus one are used
#'  
#' @details 
#' Different quality criteria used to select the best models among a list of 
#' models with different interval definitions. Quality criteria are: 
#' log-likelihood (\code{loglik}), cross-validation error as provided by the
#' function \code{\link[glmnet]{glmnet}}, two versions of the AIC (\code{AIC} 
#' and \code{AIC2}) and of the BIC (\code{BIC} and \code{BIC2}) in which the 
#' number of parameters is either the number of non null intervals or the 
#' number of non null parameters with respect to the original variables
#' 
#' @author Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}
#' 
#' @references Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016) 
#' Interpretable sparse SIR for digitized functional data.
#' \emph{Statistics and Computing}, \strong{29}(2), 255--267.
#' 
#' @seealso \code{\link{ridgeSIR}}, \code{\link{sparseSIR}}
#' 
#' @examples
#' set.seed(1140)
#' tsteps <- seq(0, 1, length = 200)
#' nsim <- 100
#' simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
#' x <- t(replicate(nsim, simulate_bm()))
#' beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
#' beta[((tsteps < 0.2) | (tsteps > 0.5)), 1] <- 0
#' beta[((tsteps < 0.6) | (tsteps > 0.75)), 2] <- 0
#' y <- log(abs(x %*% beta[ ,1]) + 1) + sqrt(abs(x %*% beta[ ,2]))
#' y <- y + rnorm(nsim, sd = 0.1)
#' res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = 10^8)
#' \donttest{res_fused <- SISIR(res_ridge, rep(1, ncol(x)), ncores = 2)
#' res_fused}
#' 
#' @return S3 object of class \code{SISIR}: a list consisting of
#' \describe{
#'    \item{\code{sEDR}}{ the estimated EDR spaces (a list of p x d matrices)}
#'    \item{\code{alpha}}{ the estimated shrinkage coefficients (a list of 
#'    vectors)}
#'    \item{\code{intervals}}{ the interval lengths (a list of vectors)}
#'    \item{\code{quality}}{ a data frame with various qualities for the model.
#'    The chosen quality measures are the same than for the function 
#'    \code{\link{sparseSIR}} plus the number of intervals \code{nbint}}
#'    \item{\code{init_sel_prop}}{ initial fraction of the coefficients which 
#'    are considered as strong zeros or strong non zeros}
#'    \item{\code{rSIR}}{ same as the input \code{object}}
#'  }
#'  
SISIR <- function(object, inter_len = rep(1, nrow(object$EDR)), sel_prop = 0.05,
                  itermax = Inf, minint = 2, parallel = TRUE, ncores = NULL) {
  
  init_sel_prop <- sel_prop
  
  nbiter <- 1
  loglik <- AIC <- BIC <- AIC2 <- BIC2 <- CVerror <- NULL
  sEDR <- list()
  alpha <- list()
  intervals <- list()
  intervals[[1]] <- inter_len
  while ((sel_prop < 0.49) & (nbiter < itermax) & length(inter_len) > minint) {
    sparse_res <- sparseSIR(object, inter_len, adaptive = TRUE, sel_prop,
                            parallel, ncores)
    
    # collect results
    CVerror <- c(CVerror, sparse_res$quality["CVerror"])
    AIC <- c(AIC, sparse_res$quality["AIC"])
    BIC <- c(BIC, sparse_res$quality["BIC"])
    AIC2 <- c(AIC2, sparse_res$quality["AIC2"])
    BIC2 <- c(BIC2, sparse_res$quality["BIC2"])
    loglik <- c(loglik, sparse_res$quality["loglik"])

    sEDR[[nbiter]] <- sparse_res$sEDR
    alpha[[nbiter]] <- sparse_res$alpha
    
    ## merge intervals (sequential from left to right)
    liste_strong0 <- which(sparse_res$adapt_res$nonzeros != 0)
    liste_strongnot0 <- which(sparse_res$adapt_res$zeros == 0)
    df_merge <- data.frame("int" = c(liste_strong0, liste_strongnot0),
                           "type" = c(rep("yes", length(liste_strong0)),
                                      rep("no", length(liste_strongnot0))))
    df_merge <- df_merge[order(df_merge$int), ]
    nb_tests <- nrow(df_merge)
    
    # merging at left border
    if (df_merge$int[1] > 1) {
      leftbd <- sum(inter_len[1:(df_merge$int[1] - 1)])
      if (leftbd < inter_len[df_merge$int[1]]) {
        inter_len[df_merge$int[1]] <- inter_len[df_merge$int[1]] + leftbd
        inter_len[1:(df_merge$int[1] - 1)] <- 0
      }
    }
    
    # merge between intervals of same types
    if (nb_tests > 1) {
      may_merge <- df_merge$type[1:(nb_tests-1)] == df_merge$type[2:nb_tests]
    } else may_merge <- FALSE
    if (sum(may_merge) > 0) {
      # merge if two intervals are consecutive
      merge1 <- (df_merge$int[1:(nb_tests-1)] + 1) == df_merge$int[2:nb_tests]
      # merge if around length is larger than between length
      around_len <- inter_len[df_merge$int[1:(nb_tests-1)]] +
        inter_len[df_merge$int[2:nb_tests]]
      between_len <- sapply(1:(nb_tests-1), function(pos) {
        sum(inter_len[(df_merge$int[pos] + 1):(df_merge$int[pos+1] - 1)])
      })
      merge2 <- around_len > between_len
      # conclusion and merge
      merge <- which(may_merge & (merge1 | merge2))
      
      # sequential merge (from left to right)
      for (pos in merge) {
        m_length <- sum(inter_len[df_merge$int[pos]:df_merge$int[pos + 1]])
        inter_len[df_merge$int[pos + 1]] <- m_length
       inter_len[df_merge$int[pos]:(df_merge$int[pos + 1] - 1)] <- 0
      }
    }
    
    # right border
    if (df_merge$int[nb_tests] < length(inter_len)) {
      rightbd <- sum(inter_len[(df_merge$int[nb_tests] + 1):length(inter_len)])
      if (rightbd < inter_len[df_merge$int[nb_tests]]) {
        inter_len[df_merge$int[nb_tests]] <- inter_len[df_merge$int[nb_tests]] +
          rightbd
        inter_len[(df_merge$int[nb_tests] + 1):length(inter_len)] <- 0
      }
    }
    
    # final intervals
    inter_len <- inter_len[inter_len != 0]
    intervals[[nbiter + 1]] <- inter_len
    cat("Current number of intervals...", length(inter_len), "\n")
    
    if (length(intervals[[nbiter + 1]]) == length(intervals[[nbiter]])) 
      sel_prop <- sel_prop + 0.01

    nbiter <- nbiter + 1
  }
  
  intervals[[nbiter]] <- NULL
  quality <- data.frame("loglik" = loglik, "AIC" = AIC, "AIC2" = AIC2, 
                        "BIC" = BIC, "BIC2" = BIC2, "CVerror" = CVerror,
                        "nbint" = unlist(lapply(intervals, length)))
  res <- list("sEDR" = sEDR, "alpha" = alpha, "intervals" = intervals,
              "quality" = quality, "init_sel_prop" = init_sel_prop, 
              "rSIR" = object)
  class(res) <- "SISIRres"
  return(res)
}

################################################################################
# Methods for objects of class SISIRres
################################################################################
#' @title Print SISIRres object
#' @name SISIRres
#' @export
#' @aliases summary.SISIRres
#' @aliases print.SISIRres
#' @description Print a summary of the result of \code{\link{SISIRres}} (
#' \code{SISIRres} object)
#' @param object a \code{SISIRres} object
#' @param x a \code{SISIRres} object
#' @param ... not used
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' @seealso \code{\link{SISIR}}
#' 
summary.SISIRres <- function(object, ...) {
  cat("Interval Sparse SIR results with:\n\n",
      object$rSIR$parameters$H, "slices\n",
      "dimension of the EDR space is:", object$rSIR$parameters$d, "\n",
      "regularization parameter is:", object$rSIR$parameters$mu2, "\n\n")
  cat("Number of fitted models:", length(object$sEDR), "\n\n")
}

#' @export
#' @rdname SISIRres
print.SISIRres <- function(x, ...) {
  summary.SISIRres(x)
}