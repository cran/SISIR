library(SISIR)

set.seed(1140)

tsteps <- seq(0, 1, length = 200)
nsim <- 100
simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
x <- t(replicate(nsim, simulate_bm()))
beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
beta[((tsteps < 0.2) || (tsteps > 0.5)), 1] <- 0
beta[((tsteps < 0.6) || (tsteps > 0.75)), 2] <- 0
y <- log(abs(x %*% beta[ ,1]) + 1) + sqrt(abs(x %*% beta[ ,2])) + 
  rnorm(nsim, sd = 0.1)

# perform ridge SIR
res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = 10^5)
# perform sparse SIR
res_sparse <- sparseSIR(res_ridge, rep(1, ncol(x)))
print(res_sparse)
# interval fusion
res_fuse <- SISIR(res_ridge, rep(1, ncol(x)), parallel = TRUE)
print(res_fuse)
ind_best <- which.min(res_fuse$quality$CVerror)
cat("Best model for CV error has", res_fuse$quality$nbint[ind_best], 
    "intervals.", "\n\n")
res_best_fuse <- sparseSIR(res_ridge, res_fuse$intervals[[ind_best]])
print(res_best_fuse)
