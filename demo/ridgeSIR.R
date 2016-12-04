library(SISIR)

set.seed(1115)
tsteps <- seq(0, 1, length = 200)
nsim <- 100
simulate_bm <- function() return(c(0, cumsum(rnorm(length(tsteps)-1, sd=1))))
x <- t(replicate(nsim, simulate_bm()))
beta <- cbind(sin(tsteps*3*pi/2), sin(tsteps*5*pi/2))
y <- log(abs(x %*% beta[ ,1])) + sqrt(abs(x %*% beta[ ,2])) +
  rnorm(nsim, sd = 0.1)

# Tuning ----------------------------------------------------------------------
list_mu2 <- 10^(0:10)
listH <- c(5, 10)
list_d <- 1:4
set.seed(1129)
res_tune <- tune.ridgeSIR(x, y, listH, list_mu2, list_d, nfolds = 10, 
                          parallel = TRUE)

## for H=10
printable_res <- matrix(res_tune[res_tune$H == 10, "cverror"],
                        ncol=length(list_d), byrow = TRUE)
rownames(printable_res) <- log10(list_mu2)
colnames(printable_res) <- list_d
plot(log10(list_mu2), printable_res[ ,ncol(printable_res)], main = "",
     ylab = "CV error", xlab = expression(log[10](mu[2])), type = "b",
     pch = "+")
chosen_mu <- which.min(printable_res[ ,ncol(printable_res)])

printable_res2 <- matrix(res_tune[res_tune$H == 10, "Rd"],
                         ncol=length(list_d), byrow = TRUE)
rownames(printable_res2) <- log10(list_mu2)
colnames(printable_res2) <- list_d
plot(list_d, printable_res2[chosen_mu, ], type="b", pch="+", xlab = "d",
     ylab = expression(R(d)))
d <- 2

plot(log10(list_mu2), printable_res[ ,d], main = "",
     ylab = "CV error", xlab = expression(log[10](mu[2])), type = "b",
     pch = "+")
which.min(printable_res[ ,ncol(printable_res)]) # consistent

# Training --------------------------------------------------------------------
res_ridge <- ridgeSIR(x, y, H = 10, d = 2, mu2 = list_mu2[chosen_mu])
# > res_ridge
# Ridge SIR results with:
#   
#   10 slices
# dimension of the EDR space is: 2 
# regularization parameter is: 1e+08 
# The EDR space is in '$EDR'. First 10 rows are:
#   
#   10 x 2 Matrix of class "dgeMatrix"
# [,1]          [,2]
# [1,] 2.494088e-22  2.240430e-21
# [2,] 1.393174e-08  6.804754e-07
# [3,] 4.201512e-07  6.958648e-07
# [4,] 1.024709e-06  9.157455e-07
# [5,] 1.238248e-06  7.735075e-07
# [6,] 1.014715e-06  5.739618e-07
# [7,] 1.558069e-06 -6.855240e-07
# [8,] 1.752981e-06 -1.081445e-06
# [9,] 1.969968e-06 -2.109960e-06
# [10,] 2.580175e-06 -9.432574e-07
