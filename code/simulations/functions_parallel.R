##library(devtools)
library(parallel)
library(MCSim)
##library(ade4)
##library(vegan)
library(lme4)
##library(bbmle)
library(optimx)
##library(dplyr)
##library(MASS)


################################################################################
## Functions to fit the competing models to abundance data
################################################################################
## 1. Hypothesis of niche and neutral dynamics
m.full <- function(ab, trait, grad, site, region, spp, ...){
    niche.neutral <- glmer(ab ~ trait + grad + I(grad^2) + 
                               trait:grad + trait:I(grad^2) + 
                               (1|spp) +
                               (1|spp:region) + (1|spp:site) + (1+grad|site), ...)
    return(niche.neutral)#, 
}

## 2. Hypothesis of niche and neutral dynamics without species traits
m.full2 <- function(ab, grad, site, region, spp, ...){
    env.neutral <- glmer(ab ~ grad + I(grad^2) + 
                             (1|spp:region) + (1|spp:site) +
                             (1+grad|spp) + (1+grad|site), ...)
    return(env.neutral)
}

## 3. Solely niche dynamics
m.niche <- function(ab, trait, grad, site, spp, ...){
    niche <- glmer(ab ~ trait + grad + I(grad^2) + 
                       trait:grad + trait:I(grad^2) + 
                       (1|spp) + (1+grad|site), ...)
    return(niche)#, 
}


## 4. Solely niche dynamics without species trait (including term to catch wrong traits)
m.env <- function(ab, grad, site, spp, ...){
    env <- glmer(ab ~ grad + I(grad^2) +
                     (1+grad|spp) + (1+grad|site), ...)
    return(env) 
}

## 5. Solely neutral dynamics
m.neutral <- function(ab, site, region, spp, ...){
    neutral <- glmer(ab ~ 1 +
                        # (1|spp) 
                     + (1|site) +
                    (1|spp:region) + (1|spp:site), ...)
    return(neutral)#, 
}


## 6. Null hypothesis
# Null model in which species and sites are independent random intercepts
m.null <- function(ab, site, region, spp, ...){
    null <- glmer(ab ~ 1 + (1|region) + (1|spp) + (1|site), ...)
    return(null)
}

################################################################################
## Functions to generate data and run the simulations
################################################################################

#' @param Nsites number of sites in the simulated metaccomunities. Should be a multiple of Nregions.
#' @param Nregions number of regions to aggregate the sites.
#' @param Nspp number of species in the simulated metacommunity
#' @param JM total number of individuals in the metacommunity
#' @param m migration parameter, to be assed to MCSim::fn.make.landscape
#' @param S.eff sampling effort, the proportion of total number of individuals in the matacoomunity that is expected to be sampled.
#' @param Nrep number of replicates of simulations of the metacommunity. The abundances in the simulated metacommunity will be the mean of the abundances over replicates.
generate.data <- function(Nsites=30, Nregions=3, Nspp=153, JM=1e6, m=0.5, Nrep=10, S.eff=0.02){
    ## Sites attributes ##
    ## Here are the xy-coordinates of sites: 3 regions with 10 sites each
    ## Distances between regions is an order of magnitude distance bewteen sites within regions
    sites <- expand.grid(x = seq(0,120, length=Nregions), y=seq(1,5, length=Nsites/Nregions))
    sites$x <- jitter(sites$x)
    sites$y <- jitter(sites$y)
    ## Each set of 10 points at same x coordinate is labelled as a region
    sites$region <- rep(letters[1:Nregions], Nsites/Nregions)
    ## Enviromental variable: 10 states, that repeat at each site (e.g. altitude)
    sites$env <- rep(seq(1,5, length=Nsites/Nregions), each=Nregions)
    ## Calculate niches for species: optimal values along the enviromental variable
    sp.opt <- runif(Nspp, min = 1, max = 5)
    ## Initial condition ##
    ## Initial condition: matrix of sites x spp
    ## Random values of species abundances 
    m0b <- matrix(rlnorm(Nspp*Nsites), Nsites, Nspp)
    ## Round values of species abundance to represent discrete values of number of individuals
    m0b <- round(m0b)
    ## Splitting species in 3 fractions that are exclusive of each region
    R.ind <- sample(letters[1:Nregions], Nspp, replace=TRUE)
    for(i in letters[1:Nregions])
        m0b[sites$region==i,R.ind!=i] <- 0
    ## Calculating Relative abundances
    m0b <- sweep(m0b, 1, apply(m0b,1,sum), FUN="/")
    ## We set m=0.5 to allow some dispersal limitation.
    ## Still a half of the deaths are replaced by locals
    simulation_landscape_det <- MCSim::fn.make.landscape(
                                           site.coords = sites[,1:2],
                                           Ef = sites$env,
                                           m = m, 
                                           JM = JM)
    ## Deteministic community ##
    ## Data frames to store simulations
    id <- data.frame(site=rep(1:Nsites, Nspp),
                     sites[,3:4],
                     spp=rep(paste("sp",1:Nspp, sep="."), each=Nsites))
    ## To store simulation results
    Nrep <- 10
    det.resu <- matrix(NA, nrow=nrow(id), ncol=Nrep)
    
    for(i in 1:Nrep){
        simu.det <- MCSim::fn.metaSIM( # simulation of deterministic community w correctly observed traits
                               landscape = simulation_landscape_det,
                               trait.Ef = sp.opt,
                               trait.Ef.sd = 0.5,
                               J.t0 = m0b,
                               n.timestep = 100, # increased to 100 time steps; initial conditions seems to persist after t=30
                               W.r = 0,
                               nu = 0,
                               speciation.limit = 0,
                               save.sim = FALSE
                           )
        det <- subset(simu.det$J.long, timestep=="100")[,4]
        det.resu[,i] <- det 
    }
    ## Abundances of each species as the mean of the 10 runs
    det <- data.frame(id, count=rowMeans(det.resu))
    ## Stochastic community ##
    simulation_landscape_sto <- MCSim::fn.make.landscape(
                                           site.coords = sites[,1:2],
                                           Ef = sites$env,
                                           m = 0.5, 
                                           JM = JM)
    sto.resu <- matrix(NA, nrow=nrow(id), ncol=Nrep)    
    for(i in 1:Nrep){
        simu.sto <- MCSim::fn.metaSIM(
                               landscape = simulation_landscape_sto,
                               trait.Ef = sp.opt,
                               trait.Ef.sd = 1000, # niche deviation changed for neutral dynamics
                               J.t0 = m0b,
                               n.timestep = 100, 
                               W.r = 200, # Dispersal Kernel no longer flat
                               nu = 0,
                               speciation.limit = 0,
                               save.sim = FALSE
                           )
        sto <- subset(simu.sto$J.long, timestep=="100")[,4]
        sto.resu[,i] <- sto 
    }
    sto <- data.frame(id, count=rowMeans(sto.resu))
    ## Poisson samples of a fraction S.eff of the total number of individuals##
    ## Deterministic community
    det.pois <- rpois(length(det$count), det$count*S.eff)
    ## Stochastic community
    sto.pois <- rpois(length(sto$count), sto$count*S.eff)
    ## Binding all data togheter ##
    data <- data.frame(site=id[,"site"], spp=id[,"spp"],
                       det=det.pois, sto=sto.pois)                   
    ## Adding species traits, gradient and spacial info
    ## Traits 
    ## A vector with wrong traits with correlation of less than 0.01 with the true traits
    cor.t <- 1
    while(cor.t>0.01){
        wrong.t <- runif(length(sp.opt), min(sp.opt), max(sp.opt))
        cor.t <- abs(cor(wrong.t, sp.opt))
    }
    ## Trait data
    trait.data <- data.frame(spp=unique(id$spp), trait= scale(sp.opt),
                         trait.wr=scale(wrong.t)) # creating vector w/ wrong traits
    ## Gradient
    env.data <- data.frame(site=unique(id$site), grad=scale(sites$env),
                           region=sites$region)
    
    ## Preparing data table for model selection
    all.data <- merge(data, env.data, by=c("site"))
    all.data <- merge(all.data, trait.data, by=c("spp"))
    return(all.data)
}

################################################################################
## Functions to fit all competing models to data generated by the function above
################################################################################
#' @param data object with simulated data
#' @param ab character, name of the variable in data that has the species abundances ("det" for abundances for deterministic simulation and "sto" for stochastic simulation)
#' @param trait character, name of the variable in data that has the species traits ("trait" for traits correlated with the abundances of species and "trait.wr" for traits uncorrelated, that is uninformative traits).
fit.them.all <- function(objeto, ab, trait){
    ## With traits
    nineu.t <- m.full(ab=objeto[,ab],
                           trait=objeto[,trait],
                           grad=objeto$grad,
                           site=objeto$site,
                           region=objeto$region, 
                           spp=objeto$spp,
                           family="poisson",
                           control=glmerControl(optimizer="bobyqa", 
                                                optCtrl=list(maxfun=5e7))
                           )
    ## Without traits
    envneu <- m.full2(ab=objeto[,ab],
                          grad=objeto$grad,
                          site=objeto$site,
                          region=objeto$region, 
                          spp=objeto$spp,
                          family="poisson",
                          control=glmerControl(optimizer="bobyqa", 
                                               optCtrl=list(maxfun=5e7))
                          )
    ## Niche dynamics
    ## With traits
    niche.t <- m.niche(ab=objeto[,ab],
                            trait=objeto[,trait],
                            grad=objeto$grad,
                            site=objeto$site, 
                            spp=objeto$spp,
                            family="poisson",
                            control=glmerControl(optimizer="bobyqa", 
                                                 optCtrl=list(maxfun=5e7))
                            )
    ## W/o traits
    env <- m.env(ab=objeto[,ab],
                     grad=objeto$grad,
                     site=objeto$site,
                     spp=objeto$spp,
                     family="poisson",
                     control=glmerControl(optimizer="bobyqa", 
                                          optCtrl=list(maxfun=5e7))
                     )

    ## Neutral dynamics
    neu <- m.neutral(ab=objeto[,ab],
                         site=objeto$site,
                         region=objeto$region, 
                         spp=objeto$spp,
                         family="poisson",
                         control=glmerControl(optimizer="bobyqa", 
                                              optCtrl=list(maxfun=5e7))
                         )

    ## Null hypothesis
    null <- m.null(ab=objeto[,ab], 
                       site=objeto$site,
                       region=objeto$region, 
                       spp=objeto$spp,
                       family="poisson",
                       control=glmerControl(optimizer="bobyqa", 
                                            optCtrl=list(maxfun=5e7))
                       )
    list(nineu.t=nineu.t, envneu=envneu, niche.t=niche.t, env=env, neu=neu, null=null)
}

################################################################################
## Parallel version of replicate
## https://rdrr.io/github/grayclhn/dbframe-R-library/man/RepParallel.html
RepParallel <- function(n, expr, simplify = "array",...) {
    answer <-
        mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
    if (!identical(simplify, FALSE) && length(answer)) 
        return(simplify2array(answer, higher = (simplify == "array")))
    else return(answer)
}
