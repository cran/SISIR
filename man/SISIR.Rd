% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseSIR.R
\name{SISIR}
\alias{SISIR}
\title{Interval Sparse SIR}
\usage{
SISIR(
  object,
  inter_len = rep(1, nrow(object$EDR)),
  sel_prop = 0.05,
  itermax = Inf,
  minint = 2,
  parallel = TRUE,
  ncores = NULL
)
}
\arguments{
\item{object}{an object of class \code{ridgeRes} as obtained from the 
function \code{\link{ridgeSIR}}}

\item{inter_len}{(numeric) vector with interval lengths for the initial 
state. Default is to set one interval for each variable (all intervals have
length 1)}

\item{sel_prop}{fraction of the coefficients that will be considered as 
strong zeros and strong non zeros. Default to 0.05}

\item{itermax}{maximum number of iterations. Default to Inf}

\item{minint}{minimum number of intervals. Default to 2}

\item{parallel}{whether the computation should be performed in parallel or
not. Logical. Default is FALSE}

\item{ncores}{number of cores to use if \code{parallel = TRUE}. If left to 
NULL, all available cores minus one are used}
}
\value{
S3 object of class \code{SISIR}: a list consisting of
\describe{
   \item{\code{sEDR}}{ the estimated EDR spaces (a list of p x d matrices)}
   \item{\code{alpha}}{ the estimated shrinkage coefficients (a list of 
   vectors)}
   \item{\code{intervals}}{ the interval lengths (a list of vectors)}
   \item{\code{quality}}{ a data frame with various qualities for the model.
   The chosen quality measures are the same than for the function 
   \code{\link{sparseSIR}} plus the number of intervals \code{nbint}}
   \item{\code{init_sel_prop}}{ initial fraction of the coefficients which 
   are considered as strong zeros or strong non zeros}
   \item{\code{rSIR}}{ same as the input \code{object}}
 }
}
\description{
\code{SISIR} performs an automatic search of relevant intervals
}
\details{
Different quality criteria used to select the best models among a list of 
models with different interval definitions. Quality criteria are: 
log-likelihood (\code{loglik}), cross-validation error as provided by the
function \code{\link[glmnet]{glmnet}}, two versions of the AIC (\code{AIC} 
and \code{AIC2}) and of the BIC (\code{BIC} and \code{BIC2}) in which the 
number of parameters is either the number of non null intervals or the 
number of non null parameters with respect to the original variables
}
\examples{
set.seed(1140)
tsteps <- seq(0, 1, length = 200)
nsim <- 100
simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
x <- t(replicate(nsim, simulate_bm()))
beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
beta[((tsteps < 0.2) | (tsteps > 0.5)), 1] <- 0
beta[((tsteps < 0.6) | (tsteps > 0.75)), 2] <- 0
y <- log(abs(x \%*\% beta[ ,1]) + 1) + sqrt(abs(x \%*\% beta[ ,2]))
y <- y + rnorm(nsim, sd = 0.1)
res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = 10^8)
\donttest{res_fused <- SISIR(res_ridge, rep(1, ncol(x)), ncores = 2)
res_fused}

}
\references{
Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016) 
Interpretable sparse SIR for digitized functional data.
\emph{Statistics and Computing}, \strong{29}(2), 255--267.
}
\seealso{
\code{\link{ridgeSIR}}, \code{\link{sparseSIR}}
}
\author{
Victor Picheny, \email{victor.picheny@inrae.fr}\cr
Remi Servien, \email{remi.servien@inrae.fr}\cr
Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}
}
