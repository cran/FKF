## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
                      collapse = TRUE,
                      comment = "#>"
)

## ---- load_package------------------------------------------------------------
library('FKF')

## ----ex1_setup----------------------------------------------------------------
n <- 1000

## Set the AR parameters
ar1 <- 0.6
ar2 <- 0.2
ma1 <- -0.2
sigma <- sqrt(0.2)

a <- arima.sim(model = list(ar = c(ar1, ar2), ma = ma1), n = n,
               innov = rnorm(n) * sigma)

## ----ex1_ss-------------------------------------------------------------------
arma21ss <- function(ar1, ar2, ma1, sigma) {
    Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
    Zt <- matrix(c(1, 0), ncol = 2)
    ct <- matrix(0)
    dt <- matrix(0, nrow = 2)
    GGt <- matrix(0)
    H <- matrix(c(1, ma1), nrow = 2) * sigma
    HHt <- H %*% t(H)
    a0 <- c(0, 0)
    P0 <- matrix(1e6, nrow = 2, ncol = 2)
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
}

## ----ex1_obj------------------------------------------------------------------
objective <- function(theta, yt) {
    sp <- arma21ss(theta["ar1"], theta["ar2"], theta["ma1"], theta["sigma"])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
    return(-ans$logLik)
}

## ----ex1_opt------------------------------------------------------------------
theta <- c(ar = c(0, 0), ma1 = 0, sigma = 1)
fit <- optim(theta, objective, yt = rbind(a), hessian = TRUE)

## ----ex1_ci-------------------------------------------------------------------
## Confidence intervals
p <- cbind(actual = c(ar1=ar1,ar2=ar2,ma1=ma1,sigma=sigma),
           estimate = fit$par,
           lowerCI = fit$par - qnorm(0.975) * sqrt(diag(solve(fit$hessian))),
           upperCI = fit$par + qnorm(0.975) * sqrt(diag(solve(fit$hessian))))
p

## ----filter-------------------------------------------------------------------
sp <- arma21ss(fit$par["ar1"], fit$par["ar2"], fit$par["ma1"], fit$par["sigma"])
ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(a))

## ----pred---------------------------------------------------------------------
plot(ans, at.idx = 1, att.idx = NA, CI = NA)
lines(a, lty = "dotted")

## ----filt---------------------------------------------------------------------
plot(ans, at.idx = NA, att.idx = 1, CI = NA)
lines(a, lty = "dotted")

## ----guassian-----------------------------------------------------------------
plot(ans, type = "resid.qq")

## ----acf----------------------------------------------------------------------
plot(ans, type = "acf")

## ----fks----------------------------------------------------------------------
sm <- fks(ans)

## ----fks_plot-----------------------------------------------------------------
plot(sm)
lines(a,col="black")

