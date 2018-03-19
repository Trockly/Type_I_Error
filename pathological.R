source("sim.R")

# num of sim. 
Nsim <- 10000
Nsubj <- 50
Nitem <- 20

data <- mkdata(Nsubj, Nitem)

#
# Type-I error: there is no effect, what is the Type-I error-rate?
# Means: b_0 = b_1 = 2000
data$y <- simdat(data, c(2000,0), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6)
models <- fitModels(data)
res1 <- lapply(1:Nsim, function(n) simstep(models, data, c(2000,0), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6))
save(res1, file="sims_path_H0_N50_M20_RII100.rda")

#
# Type-II error: there is an effect, what is the type-II error-rate 
# Means: b_0 = 2000; b_1 = 2100
data$y <- simdat(data, c(2000,25), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6)
models <- fitModels(data)
# Type-II error: there is an effect, what is the type-II error-rate 
res2 <- lapply(1:Nsim, function(n) simstep(models, data, c(2000,25), sig=300, sigRSI=100, sigRSS=0, sigRII=100, sigRIS=0, rho=0.6))
save(res2, file="sims_path_H1_N50_M20_C25_RII100.rda")
