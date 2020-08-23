# Running models and analysis for all ------------------------------------------

#
library(parallel)
source("code/simulations/functions_parallel.R")
mc.cores = detectCores()

# Already loading the simulated data -------------------------------------------
# simulated.data <- RepParallel(n=100,
#                               generate.data(),
#                               mc.cores=mc.cores, simplify=FALSE)
# save(simulated.data, file="simulated_data.RData")
load("results/RData/simulated_data.RData")

# Making models ----------------------------------------------------------------
# here only for det_rt, the only one that RData is not good
source("code/simulations/simulations.R")

# Making analysis --------------------------------------------------------------
source("code/simulations/simulation_analyses.R")

