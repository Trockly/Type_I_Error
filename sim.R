# Just creates table of covariates
mkdata <- function(N, M, C=c(-0.5,0.5)) {
  expand.grid(C=C, I=factor(1:M, labels=letters[1:M]),
              S=factor(1:N), KEEP.OUT.ATTRS = FALSE)
}

# Sampler
simdat <- function(data, bb, sig, sigRSI, sigRSS, sigRII, sigRIS, rho, nsim) {
  stopifnot(is.numeric(C <- data$C),
            is.factor(S <- data$S),
            is.factor(I <- data$I))
  # Assemble (lower) Cholesky factor of RE cov. matrices
  beta  <- c(bb[[1]], bb[[2]]); 
  names(beta) <- c("(Intercept)", "C");
  theta <- c(sigRSI/sig, sigRSS*rho/sig, sigRSS*sqrt(1-rho^2)/sig,
             sigRII/sig, sigRIS*rho/sig, sigRIS*sqrt(1-rho^2)/sig);
  names(theta) <- c("S.(Intercept)", "S.C.(Intercept)", "S.C",
                    "I.(Intercept)", "I.C.(Intercept)", "I.C");
  # Create array of simulated data
  simulate(y ~ 1 + C + (1+C|S) + (1+C|I), 
           newdata=data, family=gaussian,
           newparams=list(beta=beta, theta=theta, sigma=sig), nsim=nsim)
}

# Fit all models
fitModels <- function(data) {
  list(
    mods=list(
      m1 = lmer(y ~ 1+C + (1+C|S) + (1+C|I), data, REML=FALSE),
      m2 = lmer(y ~ 1+C + (1+C||S) + (1+C||I), data, REML=FALSE),
      m3 = lmer(y ~ 1+C + (1+C||S) + (1|I), data, REML=FALSE),
      m4 = lmer(y ~ 1+C + (1|S) + (1+C||I), data, REML=FALSE),
      m5 = lmer(y ~ 1+C + (1|S) + (1|I), data, REML=FALSE)),
    mods0=list(
      m1 = lmer(y ~ 1 + (1+C|S) + (1+C|I), data, REML=FALSE),
      m2 = lmer(y ~ 1 + (1+C||S) + (1+C||I), data, REML=FALSE),
      m3 = lmer(y ~ 1 + (1+C||S) + (1|I), data, REML=FALSE),
      m4 = lmer(y ~ 1 + (1|S) + (1+C||I), data, REML=FALSE),
      m5 = lmer(y ~ 1 + (1|S) + (1|I), data, REML=FALSE)) )  
}

# Update models given reponse
updateModels <- function(models, response) {
  list(
    mods=list(
      m1 = refit(models$mods$m1, response),
      m2 = refit(models$mods$m2, response),
      m3 = refit(models$mods$m3, response),
      m4 = refit(models$mods$m4, response),
      m5 = refit(models$mods$m5, response)),
  mods0=list(
      m1 = refit(models$mods0$m1, response),
      m2 = refit(models$mods0$m2, response),
      m3 = refit(models$mods0$m3, response),
      m4 = refit(models$mods0$m4, response),
      m5 = refit(models$mods0$m5, response)) )
}


# get t-value of the "slope" fixed effect
tval <- function(fm) coef(summary(fm))["C","t value"]
# get estimated value for the "slope"
mval <- function(fm) coef(summary(fm))["C","Estimate"]
# get std. err. for the est. value of the "slope"
sval <- function(fm) coef(summary(fm))["C", "Std. Error"]
# get smallest (relative) singular value of the RE cov.
smin <- function(fm) min(vapply(summary(rePCA(fm)), function(comp) min(comp$importance[2,]),1))


# get model stats
simstats <- function(models) {
  cbind(tstats = vapply(models$mods,tval,1),
        means  = vapply(models$mods,mval,1),
        stderr = vapply(models$mods,sval,1),
        Dev = vapply(models$mods, deviance, 1),        
        chi2 = (vapply(models$mods0, deviance, 1) - 
                  vapply(models$mods, deviance, 1)),
        smin = vapply(models$mods, smin, 1), 
        AIC = vapply(models$mods,AIC,1),
        BIC = vapply(models$mods,BIC,1))
}

# Perform a single simulation step given the "prefitted" models
simstep <- function(models, data, bb, sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6) {
  failed <- TRUE; res <- NULL;
  # A call/cc using withRestart() may be more beautiful.
  #  However, I do not know how the call/cc would affect the RNG.
  while(failed) {
    # sample responses
    resp <- simdat(data, bb, sig, sigRSI, sigRSS, sigRII, sigRIS, rho)
    tryCatch( {
      # refit models
      mods <- updateModels(models, resp);
      # get stats
      res <- simstats(mods);
      # done
      failed <- FALSE; models <- mods;
    }, 
    error=function(e) { print("err!"); failed<-TRUE; }, 
    warning=function(w) { print("oha!"); failed<-TRUE; })    
  }
  return( res );
}


barr_modsel_backward <- function(stat) {
  # Backward scheme: 
  #  Choose the more complex model only if there is 
  #  a sig. difference between the two models (LRT, alpha = 0.20)
  alpha <- 0.20
  s1 <- qchisq(1-alpha, df=1)
  s2 <- qchisq(1-alpha, df=2)
  
  # if model 1 is better than 2 -> choose 1
  if ( (stat[2,"Dev"]-stat[1,"Dev"])>s2 ) { return(1); }
  # if model 2 is better than either 3 or 4 -> choose 2
  if ( ((stat[3,"Dev"]-stat[2,"Dev"])>s1) || ((stat[4,"Dev"]-stat[2,"Dev"])>s1) ) { return(2); }
  # else if model 3 is better than 5  -> choose 3
  if ( (stat[5,"Dev"]-stat[3,"Dev"])>s1 ) { return(3); }
  # else if model 4 is better than 5  -> choose 4
  if ( (stat[5,"Dev"]-stat[4,"Dev"])>s1 ) { return(4); }
  # else choose 5
  return(5);
}


