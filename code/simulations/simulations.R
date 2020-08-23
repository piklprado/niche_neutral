# library(parallel)
# source("code/simulations/functions_parallel.R")
# mc.cores = detectCores()
#
# simulated.data <- RepParallel(n=100,
#                               generate.data(),
#                               mc.cores=mc.cores, simplify=FALSE)
# save(simulated.data, file="simulated_data.RData")
# load("results/simulated_data.RData")

## Deterministic community with correct traits
f1 <- function(x){
    result <- try(fit.them.all(x, ab = "det", trait = "trait"))
    if (class(result) == "try-error")
        return(NA)
    else
        return(result)
}

start.time <- Sys.time()
det.rt <- mclapply(simulated.data, f1, mc.cores = mc.cores)
save(det.rt, file = "simulations/RData/det_rt.RData")
rm(det.rt)
end.time <- Sys.time()
t.det.rt <- end.time - start.time
#save.image()

# Deterministic community with wrong traits
f2 <- function(x){
    result <- try(fit.them.all(x, ab = "det", trait = "trait.wr"))
    if (class(result) == "try-error")
        return(NA)
    else
        return(result)
    }
det.wt <- mclapply(simulated.data, f2, mc.cores = mc.cores)
save(det.wt, file = "simulations/RData/det_wt.RData")
rm(det.wt)
#save.image()

## Stochastic community with righ traits
f3 <- function(x){
    result <- try(fit.them.all(x, ab = "sto", trait = "trait"))
    if (class(result) == "try-error")
        return(NA)
    else
        return(result)
    }
sto.rt <- mclapply(simulated.data, f3, mc.cores = mc.cores)
save(sto.rt, file = "simulations/RData/sto_rt.RData")


write.csv(data.frame(type = c("det.rt"), duration = t.det.rt),
          "results/time_sumulations.csv")
