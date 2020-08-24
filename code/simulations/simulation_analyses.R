source("code/functions.R")
library(MuMIn)

## Deterministic community, right traits ##
## Loads the model list
load("results/RData/det_rt2.RData")
## Model selection
## Auxiliary function
f1 <- function(x, n = 6) {
    if (is.null(x))
        y <- rep(NA, n)
    else
        y <- sapply(x, MuMIn::AICc)
    y - min(y)
}
## Delta-AICc for each model in each simulated data set
det.rt.aic <- sapply(det.rt, f1)
## Proportion of the simulations in which each model was among the selected models
df.dt.rt <- apply(det.rt.aic < 2, 1, function(x) sum(x, na.rm = TRUE)/sum(!is.na(x)))

## Pseudo-R2
## Auxiliary function
f2 <- function(x) {
    if (is.null(x))
        y <- rep(NA, 7)
    else
        y <- with(x, r2.full(nineu.t, null.model = neu)$full.table[, 2])
    return(y)
}
## R2 calculations
det.rt.R2 <- sapply(det.rt, f2)
## ugly tweak to add names to the rows
rownames(det.rt.R2) <- with(det.rt[[1]], as.character(r2.full(nineu.t, null.model = neu)$full.table[,1]))
## Mean R2 for each component
df2.dt.rt <- rbind(
    apply(det.rt.R2, 1, mean),
    ## Empirical 95% CI
    apply(det.rt.R2, 1, quantile, c(0.025, 0.975))
)

write.csv(df.dt.rt, "results/simulations/dt_rt_prop.csv")
write.csv(df2.dt.rt, "results/simulations/dt_rt_effects.csv")
## To save RAM
rm(det.rt)
#save.image()


## Deterministic community, wrong traits ##
## Loads the model list
load("results/RData/det_wt2.RData")
## Model selection
## Delta-AICc for each model in each simulated data set
det.wt2.aic <- sapply(det.wt2, f1)
## Proportion of the simulations in which each model was among the selected models
df.dt.wt <- apply(det.wt2.aic < 2, 1, function(x) sum(x, na.rm = TRUE)/sum(!is.na(x)))
## Pseudo-R2
## Auxiliary function
f3 <- function(x) {
    if (is.null(x))
        y <- rep(NA, 7)
    else
        y <- with(x, r2.full2(envneu, null.model = neu)$full.table[,2])
    return(y)
}
det.wt2.R2 <- sapply(det.wt2, f3)
## ugly tweak to add names to the rows
rownames(det.wt2.R2) <- with(det.wt2[[1]], as.character(r2.full2(envneu, null.model = neu)$full.table[,1]))
## Mean R2 for each component
df2.dt.wt <- rbind(
    apply(det.wt2.R2, 1, mean, na.rm = TRUE),
    apply(det.wt2.R2, 1, quantile, c(0.025, 0.975), na.rm=TRUE)
)
## To save RAM
write.csv(df.dt.wt, "results/simulations/dt_wt_prop.csv")
write.csv(df2.dt.wt, "results/simulations/dt_wt_effects.csv")

rm(det.wt2)
#save.image()

## Stochastic community, right traits ##
## Loads the model list
load("results/RData/sto_rt.RData")
## Model selection
## Delta-AICc for each model in each simulated data set
sto.rt.aic <- sapply(sto.rt, f1)
## Proportion of the simulations in which each model was among the selected models
df.sto.rt <- apply(sto.rt.aic < 2, 1, function(x) sum(x, na.rm = TRUE)/sum(!is.na(x)))
## Pseudo-R2
## Auxiliary function
f4 <- function(x) {
    if (is.null(x))
        y <- rep(NA, 6)
    else
        y <- with(x, r2.neutral(neu, null.model = neu)$full.table[, 2])
    return(y)
}
sto.rt.R2 <- sapply(sto.rt, f4)
## ugly tweak to add names to the rows
rownames(sto.rt.R2) <- with(sto.rt[[1]], as.character(r2.neutral(neu, null.model = neu)$full.table[,1]))
## Mean R2 for each component
df2.sto.rt <- rbind(
    apply(sto.rt.R2, 1, mean, na.rm = TRUE),
    apply(sto.rt.R2, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
)
## To save RAM
write.csv(df.sto.rt, "results/simulations/sto_rt_prop.csv")
write.csv(df2.sto.rt, "results/simulations/sto_rt_effects.csv")

rm(sto.rt)
#save.image()
