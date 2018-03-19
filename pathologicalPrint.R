load("sims_path_H0_N50_M20_RII100.rda")
load("sims_path_H1_N50_M20_C25_RII100.rda")

# Get false-negative rate at about 0.05 sig. level for each case
N <- length(res1)
counts <- rowSums(sapply(res1, function(fit) { fit[,"chi2"] > qchisq(0.95, 1) } )); 
# get stats
for (i in 1:5) {
  pt <- prop.test(counts[[i]], N, p=0.05);
  print(sprintf("Type-I model %s: %s (%s, %s):", i, pt$estimate, pt$conf.int[[1]], pt$conf.int[[2]]))
}

counts <- rowSums(sapply(res1, function(fit) { fit[,"chi2"] > qchisq(0.95, 1) } )); 
# get stats
for (i in 1:5) {
  pt <- prop.test(counts[[i]], N, p=0.05);
  print(sprintf("Power model %s: %s (%s, %s):", i, pt$estimate, pt$conf.int[[1]], pt$conf.int[[2]]))
}
