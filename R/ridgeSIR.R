#' SI-SIR
#'
################################################################################
# ridge SIR
################################################################################
#' @import foreach
#' @import doParallel
#' @importFrom Matrix forceSymmetric
#' @importFrom stats cov
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom parallel detectCores
#' @importFrom RSpectra eigs
#' @importFrom expm sqrtm
#' 
#' @title ridge SIR
#' @export
#'
#' @description
#' \code{ridgeSIR} performs the first step of the method (ridge regularization
#' of SIR)
#'
#' @param x explanatory variables (numeric matrix or data frame)
#' @param y target variable (numeric vector)
#' @param H number of slices (integer)
#' @param mu2 ridge regularization parameter (numeric, positive)
#' @param d number of dimensions to be kept
#' 
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' 
#' @references {Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016) 
#' Interpretable sparse SIR for digitized functional data. 
#' \emph{Statistics and Computing}, \strong{29}(2), 255--267.}
#' 
#' @seealso \code{\link{sparseSIR}}, \code{\link{SISIR}}, 
#' \code{\link{tune.ridgeSIR}}
#' 
#' @examples
#' set.seed(1140)
#' tsteps <- seq(0, 1, length = 50)
#' simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
#' x <- t(replicate(50, simulate_bm()))
#' beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2)) 
#' y <- log(abs(x %*% beta[ ,1])) + sqrt(abs(x %*% beta[ ,2]))
#' y <- y + rnorm(50, sd = 0.1)
#' res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = 10^8)
#' \dontrun{print(res_ridge)}
#' 
#' @return S3 object of class \code{ridgeRes}: a list consisting of
#' \itemize{
#'    \item{\code{EDR}}{ the estimated EDR space (a p x d matrix)}
#'    \item{\code{condC}}{ the estimated slice projection on EDR (a d x H 
#'    matrix)}
#'    \item{\code{eigenvalues}}{ the eigenvalues obtained during the generalized 
#'    eigendecomposition performed by SIR}
#'    \item{\code{parameters}}{ a list of hyper-parameters for the method: 
#'    \itemize{
#'      \item{\code{H}}{ number of slices}
#'      \item{\code{d}}{ dimension of the EDR space}
#'      \item{\code{mu2}}{ regularization parameter for the ridge penalty}
#'    }}
#'    \item{\code{utils}}{ useful outputs for further computations:
#'    \itemize{
#'      \item{\code{Sigma}}{ covariance matrix for x}
#'      \item{\code{slices}}{ slice number for all observations}
#'      \item{\code{invsqrtS}}{ value of the inverse square root of the 
#'      regularized covariance matrix for x}
#'    }}
#'  }
#'  

ridgeSIR <- function(x, y, H, d, mu2 = NULL) {
  oldwarn <- getOption("warn")
  options(warn = -1)
  
  # precomputations
  ## for slices
  slices <- makeSlices(y, H)
  
  ## for x  
  pre <- preSIR(x, slices, H)
  
  ## for mu2
  suppressWarnings({
    smallest <- eigs(pre$cov_x, 1, which = "SR")$values
  })

  if (is.null(mu2)) {
    if (is.null(smallest)) {
      mu2 <- 10^(-5)
    } else if (smallest == Inf) {
      mu2 <- 10^(-5)
    } else if (smallest > 10^(-5)) {
      mu2 <- 0
    } else if (smallest < 0) {
      mu2 <- - smallest * 10^8
    } else mu2 <- 10^(-5)
  } else {
    if (length(smallest) > 0) {
      if ((smallest < 0) && (mu2 < 10^8 * abs(smallest))) {
        warning(paste("\n\nThe value provided for mu2 might lead to numerical",
                      "instabilities.\nMinimum advised value for mu2 is", 
                       smallest * 10^8))
      }
    }
  }
  
  # processing SIR
  resSIR <- processSIR(pre, d, mu2)
  
  # output
  precomputations <- list("Sigma" = pre$cov_x, "slices" = slices$slices,
                          "cmean_x" = pre$cmean_x, "norm_EDR" = resSIR$norm_EDR,
                          "invsqrtS" = resSIR$invsqrtS)
  orig_data <- list("x" = x, "y" = y)
  parameters <- list("H" = H, "d" = d, "mu2" = mu2)
  
  res <- list("EDR" = resSIR$EDR, "condC" = resSIR$condC,
              "eigenvalues" = resSIR$eigenvalues, "data" = orig_data,
              "parameters" = parameters, "utils" = precomputations)
  class(res) <- "ridgeRes"
  
  options(warn = oldwarn)
  return(res)
}

makeSlices <- function(y, H) {
  slices_bounds <- quantile(y, probs=seq(0, 1, length = H+1))
  slices_y <- cut(y, breaks=slices_bounds, labels=FALSE, include.lowest=TRUE)
  freq_slices <- as.vector(table(slices_y)) / length(y)
  return(list("slices" = slices_y, "freq" = freq_slices))
}

preSIR <- function(x, slices, H) {
  p <- ncol(x)
  mean_x <- colMeans(x)
  cov_x <- cov(x)
  
  ## conditional mean and variance
  cmean_x <- split(as.data.frame(x), slices$slices)
  cmean_x <- matrix(unlist(lapply(cmean_x, colMeans)), ncol = H)
  cmean_x <- cmean_x - mean_x
  ccov_x <- cmean_x %*% tcrossprod(diag(slices$freq), cmean_x)
  
  return(list("cmean_x" = cmean_x, "cov_x" = cov_x, "ccov_x" = ccov_x))
}

processSIR <- function(pre, d, mu2) {
  # computing inverse and square root of matrices
  norm_EDR <- pre$cov_x + mu2 * diag(ncol(pre$cov_x))
  suppressWarnings({sqrtS <- forceSymmetric(sqrtm(norm_EDR))})
  invsqrtS <- solve(sqrtS)
  
  # EDR
  eta <- forceSymmetric(invsqrtS %*% pre$ccov_x %*% invsqrtS)
  res <- eigs(eta, d)
  EDR <- invsqrtS %*% res$vectors
  eig_vals <- res$values
  
  # C
  cC <- crossprod(EDR, pre$cmean_x)
  
  return(list("EDR" = EDR, "condC" = cC, "eigenvalues" = eig_vals, 
              "norm_EDR" = norm_EDR, "invsqrtS" = invsqrtS))
}

################################################################################
# Methods for objects of class ridgeRes
################################################################################
#' @title Print ridgeRes object
#' @name ridgeRes
#' @export
#' @aliases summary.ridgeRes()
#' @aliases print.ridgeRes
#' @aliases ridgeRes-class
#' @description Print a summary of the result of \code{\link{ridgeSIR}} (
#' \code{ridgeRes} object)
#' @param object a \code{ridgeRes} object
#' @param x a \code{ridgeRes} object
#' @param ... not used
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr 
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' @seealso \code{\link{ridgeSIR}}
#' 
summary.ridgeRes <- function(object, ...) {
  cat("Ridge SIR results with:\n\n",
      object$parameters$H, "slices\n",
      "dimension of the EDR space is:", object$parameters$d, "\n",
      "regularization parameter is:", object$parameters$mu2, "\n")
  cat("The EDR space is in '$EDR'. First 10 rows are:\n\n")
  print(head(object$EDR, 10))
}

#' @export
#' @rdname ridgeRes
print.ridgeRes <- function(x, ...) {
  summary.ridgeRes(x)
}

################################################################################
# Cross validation for ridge SIR
################################################################################
#' @title Cross-Validation for ridge SIR
#' @export
#' @aliases tune.ridgeSIR
#'
#' @description
#' \code{tune.ridgeSIR} performs a Cross Validation for ridge SIR estimation
#'
#' @param x explatory variables (numeric matrix or data frame)
#' @param y target variable (numeric vector)
#' @param listH list of the number of slices to be tested (numeric vector)
#' @param list_mu2 list of ridge regularization parameters to be tested 
#' (numeric vector)
#' @param list_d list of the dimensions to be tested (numeric vector)
#' @param nfolds number of folds for the cross validation. Default is 10
#' @param parallel whether the computation should be performed in parallel or
#' not. Logical. Default is FALSE
#' @param ncores number of cores to use if \code{parallel = TRUE}. If left to 
#' NULL, all available cores minus one are used
#' 
#' @author {Victor Picheny, \email{victor.picheny@inrae.fr}\cr
#' Remi Servien, \email{remi.servien@inrae.fr}\cr
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#' 
#' @references {Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016) 
#' Interpretable sparse SIR for digitized functional data.
#' \emph{Statistics and Computing}, \strong{29}(2), 255--267.}
#' 
#' @seealso \code{\link{ridgeSIR}}
#' 
#' @examples
#' set.seed(1115)
#' tsteps <- seq(0, 1, length = 200)
#' nsim <- 100
#' simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
#' x <- t(replicate(nsim, simulate_bm()))
#' beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
#' y <- log(abs(x %*% beta[ ,1])) + sqrt(abs(x %*% beta[ ,2]))
#' y <- y + rnorm(nsim, sd = 0.1)
#' list_mu2 <- 10^(0:10)
#' listH <- c(5, 10)
#' list_d <- 1:4
#' set.seed(1129)
#' \dontrun{
#' res_tune <- tune.ridgeSIR(x, y, listH, list_mu2, list_d, 
#'                           nfolds = 10, parallel = TRUE)}
#' 
#' @return a data frame with tested parameters and corresponding CV error and 
#' estimation of R(d)
#' 
tune <- function(x, y, listH, list_mu2, list_d, nfolds, parallel, ncores) {
  UseMethod("tune")
}

#' @export
#' @rdname tune
tune.ridgeSIR <- function(x, y, listH, list_mu2, list_d, nfolds = 10, 
                          parallel = TRUE, ncores = NULL) {
  if (parallel) {
    if (is.null(ncores)) ncores <- min(c(detectCores() - 1, nfolds))
  }
  
  # create nfolds folds
  n <- nrow(x)
  size_fold <- round(n / nfolds)
  fold_num <- c(rep(1:(nfolds - 1), each = size_fold),
                rep(nfolds, n - (nfolds - 1) * size_fold))
  fold_num <- sample(fold_num, n, replace = FALSE)
  
  # initialization
  cv_errors <- NULL
  Rd <- NULL
  
  # loop over h
  for (h in listH) {
    print(paste("Processing h =", h))
    ## create slices
    slices <- makeSlices(y, h)
    
    ## loop over mu2 and d to obtain global projector (reference for Rd)
    max_d <- max(list_d)
    pre <- preSIR(x, slices, h)
    global_projector <- list()
    for (ind_mu2 in seq_along(list_mu2)) {
      mu2 <- list_mu2[ind_mu2]
      res_SIR <- processSIR(pre, max_d, mu2)
      global_projector[[ind_mu2]] <- sapply(list_d, function(d) {
        projectorEDR(res_SIR$EDR[ ,1:d,drop = FALSE], pre$ccov_x, mu2)
      }, simplify = FALSE)
    }
    
    ## compute this fold error
    ### parallel version
    if (parallel) {
      registerDoParallel(ncores)
      res <- foreach(cur_fold = 1:nfolds, .combine = "+",
                     .export = c("process_fold", "projectorEDR", "computeCV",
                                 "computeRd", "ridgeSIR", "preSIR",
                                 "processSIR", "makeSlices"),
                     .packages = c("expm", "RSpectra")) %dopar% {
        return(process_fold(cur_fold, x, slices$slices, fold_num, h, list_mu2,
                            list_d, global_projector))
                     }
      stopImplicitCluster()
    } else {### not parallel
      res <- foreach(cur_fold = 1:nfolds, .combine = "+") %do% {
        print(paste("process fold", cur_fold))
        return(process_fold(cur_fold, x, slices$slices, fold_num, h, list_mu2,
                            list_d, global_projector))
      }
    }
    cv_errors <- c(cv_errors, res[ ,1] / nfolds)
    Rd <- c(Rd, rep(list_d, length(list_mu2)) - res[ ,2] / nfolds)
  }
  
  lH <- length(listH)
  lmu2 <- length(list_mu2)
  ld <- length(list_d)
  out <- data.frame("H"   = rep(listH, each = lmu2 * ld), 
                    "mu2" = rep(rep(list_mu2, each = ld), lH),
                    "d"   = rep(list_d, lH * lmu2), 
                    "cverror" = cv_errors, "Rd" = Rd)
  
  return(out)
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("cur_fold"))

process_fold <- function(cur_fold, x, slices_y, fold_num, h, list_mu2, list_d,
                         global_projector) {
  # preprocessing common to all mu
  ## extract data
  test_ind <- which(fold_num == cur_fold)
  x_test <- x[test_ind, ]
  x_train <- x[- test_ind, ]
  slices_test <- slices_y[test_ind]
  freq_test <- table(slices_test) / length(slices_test)
  slices_train <- slices_y[- test_ind]
  freq_train <- table(slices_train) / length(slices_train)
  slices <- list("slices" = slices_train, "freq" = freq_train)
  ## SIR precomputations
  pre <- preSIR(x_train, slices, h)
  d_max <- max(list_d)
  
  # loop over mu2
  cv_errors <- NULL
  Rd <- NULL
  for (ind_mu2 in seq_along(list_mu2)) {
    # train on train data set
    mu2 <- list_mu2[ind_mu2]
    cur_model <- processSIR(pre = pre, d = d_max, mu2 = mu2)
    
    # compute CV error (based on MLR)
    unique_test <- sort(unique(slices_test))
    unique_train <- sort(unique(slices_train))
    mean_x <- colMeans(x_test)
    cmean_x <- lapply(split(as.data.frame(x_test), slices_test), colMeans)
    cmean_x <- matrix(unlist(cmean_x), ncol = length(unique_test))
    
    if (length(unique_test) < length(unique_train)) {
      # case where slices in test are not in train: the other case is not supported
      kept_train <- which(unique_train %in% unique_test)
    } else kept_train <- seq_along(unique_test)
    topred <- cmean_x - mean_x
    explanatory <- cov(x_test)
    smallest <- eigs(explanatory, 1, which = "SR")$value
    norm_test <- explanatory + abs(smallest) * 10^3 * diag(ncol(explanatory))
    norm_test <- forceSymmetric(norm_test)
    norm_test <- solve(norm_test)
    
    cv_errors <- c(cv_errors, sapply(list_d, function(d) {
      computeCV(topred, explanatory, freq_test, norm_test,
                cur_model$EDR[ ,1:d,drop = FALSE] %*% 
                  cur_model$condC[1:d,kept_train,drop = FALSE])
    }))
    
    # compute projection and Rd
    cur_projs <- sapply(list_d, function(d) {
      curEDR <- cur_model$EDR[ ,1:d,drop=FALSE]
      res <- projectorEDR(curEDR, pre$ccov_x, mu2)
    }, simplify = FALSE)
    Rd <- c(Rd, 
            unlist(mapply(computeRd, cur_projs, global_projector[[ind_mu2]])))
  }
  
  return(cbind(cv_errors, Rd))
}

projectorEDR <- function(EDR, ccov_x, mu2) {
  projector <- tcrossprod(EDR) %*% (ccov_x + mu2 * diag(nrow(EDR)))
  return(projector)
}

computeCV <- function(topred, explanatory, freq_test, norm, parameters) {
  quant <- topred - explanatory %*% parameters
  quant <- crossprod(quant, norm) %*% quant
  quant <- diag(quant) * freq_test
  
  return(sum(quant) / sum(freq_test))
}

computeRd <- function(projector, global) {
  sum(diag(projector %*% global))
}
