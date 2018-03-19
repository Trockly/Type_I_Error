library(mgcv)

source("sim.R")


Nsubj <- 30  # 50
Nitem <- 10  # 20 

load(sprintf("sims_scan_H0_N%i_M%i_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))
load(sprintf("sims_scan_H1_N%i_M%i_C25_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))


#
# model selection rates by AIC
#
aicmod1 <- data.frame(sig=sigs, det=sapply(res1, function(stat) (1==which.min(stat[,"AIC"]))) )
aicmod234 <- data.frame(sig=sigs, det=sapply(res1, function(stat) ((2==which.min(stat[,"AIC"])) ||
                                                                   (3==which.min(stat[,"AIC"])) ||
                                                                   (4==which.min(stat[,"AIC"])))) )
aicmod5 <- data.frame(sig=sigs, det=sapply(res1, function(stat) (5==which.min(stat[,"AIC"]))) )

aic_smooth1 <- gam(det~s(sig, k=5), family=binomial(), data=aicmod1)
aic_smooth234 <- gam(det~s(sig, k=5), family=binomial(), data=aicmod234)
aic_smooth5 <- gam(det~s(sig, k=5), family=binomial(), data=aicmod5)

# Plot model selection rates by backward LRT
bwdmod1 <- data.frame(sig=sigs, det=sapply(res1, function(stat) (1==barr_modsel_backward(stat))) )
bwdmod234 <- data.frame(sig=sigs, det=sapply(res1, function(stat) ((2==barr_modsel_backward(stat)) ||
                                                                     (3==barr_modsel_backward(stat)) ||
                                                                     (4==barr_modsel_backward(stat)))) )
bwdmod5 <- data.frame(sig=sigs, det=sapply(res1, function(stat) (5==barr_modsel_backward(stat))) )

bwd_smooth1 <- gam(det~s(sig, k=5), family=binomial(), data=bwdmod1)
bwd_smooth234 <- gam(det~s(sig, k=5), family=binomial(), data=bwdmod234)
bwd_smooth5 <- gam(det~s(sig, k=5), family=binomial(), data=bwdmod5)



#pdf(sprintf("modelsel_%i_%i_20k.pdf", Nsubj, Nitem), width=7.2, height=6)
par(mfrow=c(2,1), mar=c(.5,.5,0,0)+.1, oma=c(3,4,4,2))
plot(sigs, logistic(predict(aic_smooth1, se.fit=FALSE)), type="l", lty=1, ylim=c(0,1),
     xlab="", ylab="AIC Model selection rate", xaxt="n")
lines(sigs, logistic(predict(aic_smooth234, se.fit=FALSE)), lty=2)
lines(sigs, logistic(predict(aic_smooth5, se.fit=FALSE)), lty=3)
legend(x=0, y=.6, c("Maximal model (1)", "Reduced models (2,3,4)", "Intercept model (5)"),
       lty=c(1,2,3))
plot(sigs, logistic(predict(bwd_smooth1, se.fit=FALSE)), type="l", lty=1, ylim=c(0,1),
     xlab="", ylab="LRT backward model selection rate")
lines(sigs, logistic(predict(bwd_smooth234, se.fit=FALSE)), lty=2)
lines(sigs, logistic(predict(bwd_smooth5, se.fit=FALSE)), lty=3)
legend(x=68, y=.6, c("Maximal model (1)", "Reduced models (2,3,4)", "Intercept model (5)"),
       lty=c(1,2,3))
#dev.off()

