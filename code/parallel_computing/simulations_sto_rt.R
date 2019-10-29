library(parallel)
source("functions_parallel.R")
load("simulated_data.RData")
mc.cores=8

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

