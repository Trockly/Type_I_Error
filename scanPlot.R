library(mgcv)
library(Hmisc)

source("sim.R")

Nsubj <- 30 # 50
Nitem <- 10 # 20


load(sprintf("sims_scan_H0_N%i_M%i_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))
load(sprintf("sims_scan_H1_N30_M10_C25_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))
barr_res <- readRDS(sprintf("barr_res_%i_%i.rds", Nsubj, Nitem))

mod1res <- data.frame(
  sig=sigs, det=sapply(res1, function(stat) (stat[1,"chi2"] > qchisq(0.95, 1))) )
mod5res <- data.frame( 
  sig=sigs, det=sapply(res1, function(stat) (stat[5,"chi2"] > qchisq(0.95, 1))) )
aicres <- data.frame(
  sig=sigs, det=sapply(res1, function(stat) (stat[which.min(stat[,"AIC"]), "chi2"] > qchisq(0.95, 1))) )
backres <- data.frame(
  sig=sigs, det=sapply(res1, function(stat) (stat[barr_modsel_backward(stat), "chi2"] > qchisq(0.95, 1))) )

mod1_smooth1 <- gam(det~s(sig), family=binomial(), data=mod1res)
mod5_smooth1 <- gam(det~s(sig), family=binomial(), data=mod5res)
aic_smooth1  <- gam(det~s(sig), family=binomial(), data=aicres)
back_smooth1 <- gam(det~s(sig), family=binomial(), data=backres)

mod1_loe1     <- loess(det~sig, span=0.9, data=mod1res)
aic_loe1     <- loess(det~sig, span=0.9, data=aicres)
back_loe1     <- loess(det~sig, span=0.9, data=backres)

mod1res <- data.frame(
  sig=sigs, det=sapply(res2, function(stat) (stat[1,"chi2"] > qchisq(0.95, 1))) )
mod5res <- data.frame(
  sig=sigs, det=sapply(res2, function(stat) (stat[5,"chi2"] > qchisq(0.95, 1))) )
aicres <- data.frame(
  sig=sigs, det=sapply(res2, function(stat) (stat[which.min(stat[,"AIC"]), "chi2"] > qchisq(0.95, 1))) )
backres <- data.frame(
  sig=sigs, det=sapply(res2, function(stat) (stat[barr_modsel_backward(stat), "chi2"] > qchisq(0.95, 1))) )

mod1_smooth2 <- gam(det~s(sig), family=binomial(), data=mod1res)
mod5_smooth2 <- gam(det~s(sig), family=binomial(), data=mod5res)
aic_smooth2  <- gam(det~s(sig), family=binomial(), data=aicres)
back_smooth2 <- gam(det~s(sig), family=binomial(), data=backres)

mod1_loe2    <- loess(det~sig, span=0.9, data=mod1res)
aic_loe2     <- loess(det~sig, span=0.9, data=aicres)
back_loe2    <- loess(det~sig, span=0.9, data=backres)

logistic <- function(x) { exp(x)/(1+exp(x)) }

t1lim <- c(0.0, 0.1)
pwlim <- c(0.0, 0.5)

#pdf(sprintf("scan_compare_%i_%i_20k.pdf", Nsubj, Nitem), width=6, height=5)
par(mfrow=c(3,2), mar=c(.5,.5,0,0)+.1, oma=c(3,4,4,2))
# Type 1 error mod 1
tmp <- predict(mod1_smooth1, se.fit=TRUE)
t1max <- logistic(max(tmp$fit))
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim, xaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
barr <- barr_res[[1]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power mod 1
tmp0 <- predict(mod1_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp0$fit), type="l", lty=1, ylim=pwlim, yaxt="n", xaxt="n")
lines(sigs, logistic(tmp0$fit+2*tmp0$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit-2*tmp0$se.fit), lty=2)
barr <- barr_res[[4]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

# Type 1 error AIC selected
tmp <- predict(aic_smooth1, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim, xaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
barr <- barr_res[[2]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power AIC selected
tmp <- predict(aic_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=pwlim, xaxt="n", yaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit), lty=3)
barr <- barr_res[[5]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

# Type 1 error LRT back
tmp <- predict(back_smooth1, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim)
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
barr <- barr_res[[3]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power, LRT back
tmp <- predict(back_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=pwlim, yaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit), lty=3)
barr <- barr_res[[6]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

#dev.off()


#pdf(sprintf("scan_comploess_%i_%i_20k.pdf", Nsubj, Nitem), width=6, height=5)
par(mfrow=c(3,2), mar=c(.5,.5,0,0)+.1, oma=c(3,4,4,2))
# Type 1 error mod 1
tmp <- predict(mod1_loe1, se.fit=TRUE)
plot(sigs, tmp, type="l", lty=1, ylim=t1lim, xaxt="n")
tmp <- predict(mod1_smooth1, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
barr <- barr_res[[1]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power mod 1
tmp0 <- predict(mod1_loe2, se.fit=TRUE)
plot(sigs, tmp0, type="l", lty=1, ylim=pwlim, yaxt="n", xaxt="n")
tmp <- predict(mod1_smooth2, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
barr <- barr_res[[4]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

# Type 1 error AIC selected
tmp <- predict(aic_loe1, se.fit=TRUE)
plot(sigs, tmp, type="l", lty=1, ylim=t1lim, xaxt="n")
tmp <- predict(aic_smooth1, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
barr <- barr_res[[2]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power AIC selected
tmp <- predict(aic_loe2, se.fit=TRUE)
plot(sigs, tmp, type="l", lty=1, ylim=pwlim, xaxt="n", yaxt="n")
tmp <- predict(aic_smooth2, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
lines(sigs, tmp0, lty=3)
barr <- barr_res[[5]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

# Type 1 error LRT back
tmp <- predict(back_loe1, se.fit=TRUE)
plot(sigs, tmp, type="l", lty=1, ylim=t1lim)
tmp <- predict(back_smooth1, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
barr <- barr_res[[3]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power, LRT back
tmp <- predict(back_loe2, se.fit=TRUE)
plot(sigs, tmp, type="l", lty=1, ylim=pwlim, yaxt="n")
tmp <- predict(back_smooth2, se.fit=TRUE)
lines(sigs, logistic(tmp$fit), lty=2)
lines(sigs, tmp0, lty=3)
barr <- barr_res[[6]]
qlower <- qbeta(0.025, barr$value*10e3+1, (1-barr$value)*10e3+1)
qupper <- qbeta(0.975, barr$value*10e3+1, (1-barr$value)*10e3+1)
errbar(barr$slope, barr$value, qupper, qlower, add=T)
axis(4)

#dev.off()


#pdf(sprintf("scan_result_%i_%i_20k.pdf", Nsubj, Nitem), width=6, height=5)
par(mfrow=c(3,2), mar=c(.5,.5,0,0)+.1, oma=c(3,4,4,2))
# Type 1 error mod 1
tmp <- predict(mod1_smooth1, se.fit=TRUE)
t1max <- logistic(max(tmp$fit))
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim, xaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power mod 1
tmp0 <- predict(mod1_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp0$fit), type="l", lty=1, ylim=pwlim, yaxt="n", xaxt="n")
lines(sigs, logistic(tmp0$fit+2*tmp0$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit-2*tmp0$se.fit), lty=2)
axis(4)


# Type 1 error AIC selected
tmp <- predict(aic_smooth1, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim, xaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power AIC selected
tmp <- predict(aic_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=pwlim, xaxt="n", yaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit), lty=3)
axis(4)


# Type 1 error LRT back
tmp <- predict(back_smooth1, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=t1lim)
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
abline(h=0.05)
abline(h=t1max, lty=3)

# Power, LRT back
tmp <- predict(back_smooth2, se.fit=TRUE)
plot(sigs, logistic(tmp$fit), type="l", lty=1, ylim=pwlim, yaxt="n")
lines(sigs, logistic(tmp$fit+2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp$fit-2*tmp$se.fit), lty=2)
lines(sigs, logistic(tmp0$fit), lty=3)
axis(4)

#dev.off()
