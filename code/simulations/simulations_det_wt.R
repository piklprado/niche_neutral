library(parallel)
source("functions_parallel.R")
mc.cores=6
load("simulated_data.RData")
## Deterministic community with wrong traits
f2 <- function(x){
    result <- try(fit.them.all(x, ab="det", trait="trait.wr"))
    if(class(result)=="try-error")
        return(NA)
    else
        return(result)
    }
det.wt2 <- mclapply(simulated.data, f2, mc.cores=mc.cores)
save(det.wt2, file="det_wt2.RData")
rm(det.wt)
save.image()

