source("sim.R")

# num of sim. 
Nsim <- 20000
Nsubj <- 30  # 50
Nitem <- 10  # 20 


# Scan random slope sd. from 0 to 120
sigs <- seq(0, 120, length.out=Nsim)

data <- mkdata(Nsubj, Nitem)

# for reproducibility 
# (should not result in a sig. difference if sample size is sufficient)
set.seed(1234321)                 

# Sample a first response
data$y <- simdat(data, c(2000,0), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6)
models <- fitModels(data)
  
# Type-I error: there is no effect, what is the Type-I error-rate?
system.time(res1 <- lapply(sigs, function(sig) simstep(models, data, c(2000,0), sig=300, sigRSI=100, sigRSS=sig, sigRII=100, sigRIS=sig, rho=0.6)))
save(sigs, res1, file=sprintf("sims_scan_H0_N%i_M%i_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))

# Sample a first response
data$y <- simdat(data, c(2000,25), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6)
models <- fitModels(data)

# Type-II error: there is an effect, what is the type-II error-rate 
system.time(res2 <- lapply(sigs, function(sig) simstep(models, data, c(2000,25), sig=300, sigRSI=100, sigRSS=sig, sigRII=100, sigRIS=sig, rho=0.6)))
save(sigs, res2, file=sprintf("sims_scan_H1_N%i_M%i_C25_RSI100_RII100_rho06_20k.rda", Nsubj, Nitem))

