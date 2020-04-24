library(parallel)
source("functions_parallel.R")
mc.cores=8

simulated.data <- RepParallel(n=100,
                              generate.data(),
                              mc.cores=mc.cores, simplify=FALSE)
save(simulated.data, file="simulated_data.RData")

## Deterministic community with correct traits
f1 <- function(x){
    result <- try(fit.them.all(x, ab="det", trait="trait"))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
    }
det.rt <- mclapply(simulated.data, f1, mc.cores=mc.cores)
save(det.rt, file="det_rt.RData")
rm(det.rt)
save.image()

## Deterministic community with wrong traits
f2 <- function(x){
    result <- try(fit.them.all(x, ab="det", trait="trait.wr"))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
    }
det.wt <- mclapply(simulated.data, f2, mc.cores=mc.cores)
save(det.wt, file="det_wt.RData")
rm(det.wt)
save.image()

## Stochastic community with righ traits
f3 <- function(x){
    result <- try(fit.them.all(x, ab="sto", trait="trait"))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
    }
sto.rt <- mclapply(simulated.data, f3, mc.cores=mc.cores)
save(sto.rt, file="sto_rt.RData")

