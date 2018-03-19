# Balancing Type I error and power in linear mixed models

This archive contains the simulation code (and results) for the article

*Balancing Type I error and power in linear mixed models*
http://dx.doi.org/10.1016/j.jml.2017.01.001

## Contents
1. *sim.R* - Contains basically all simulation and model selection code needed. 
*mkdata()* - Generates a data.frame containing every combination of condition, item and subject.
*simdat()* - Samples one ore more responses given the variance/covariance parameters for every condition item and subject.
*fitModels()* - Fits all models to the data.
*updateModels()* - Re-fits all models to some new data (much faster than calling fitModels every time).
*simstats()* - Assess model-statistics for a set of fitted models.
*simstep()* - Performs a single simulation step. Samples new response, **re**fits models and obtains model statistics.
*barr_modsel_backward()* - Implements the LRT backward model selection. AIC etc. are rather trivial from the model statistics. Hence they are implemented in-line where needed.
2. *pathological.R* - Samples 10k responses from the minimal model are stores the results of all models in RDA files.
3. *scan.R* - Scans the SD of the random slopes in 20k steps and stores the results in RDA files.
4. *modselPlot.R, pathologicalPrint.R and scanPlot.R* - Output results.

## Raw simulation results
The simulations take ages, hence we added our raw simulation results to this archiv (the RDA files). I you want to alter the simulations (e.g. altering the number of items/subjects or variances) you need to re-run the simulations by calling *pathological.R or scan.R*.

