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
                         (1|spp) + (1|site) +
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
## Functions to calculate R-squared for each random effect, each function for
## Each possible best model
## adapted from Nakagawa and Johnson (2014)
## Works only for Poisson GLMMs fitted by the functions above
################################################################################


### Funtion for the full model
## PI: added 'null.model' argument and '...' argument to the fit of the null model (to allow optmization control)
r2.full <- function(model, null.model, ...){
    ## Function to calculate the null model, model with all random terms
    best.null <- function(model) {
        parens <- function(x) paste0("(",x,")")
        onlyBars <- function(form) reformulate(sapply(findbars(form),
                                                      function(x)  parens(deparse(x))),
                                               response=".")
        onlyBars(formula(model))
        best.null <- update(model,onlyBars(formula(model)), ...)
        return(best.null)
    }
    ## Calculates null model
    if(missing(null.model))
       m0 <- best.null(model)
    else
        m0 <- null.model
    ## Variance for fixed effects
    VarF <- var(as.vector(fixef(model) %*% t(model@pp$X)))
    ## Variance for slope/intercept random effects
    X <- model.matrix(model)
    n <- nrow(X)
    Z <- X[, c("(Intercept)","grad")]
    sigma <- VarCorr(model)$site
    VarSite <- sum(diag(Z %*% sigma %*% t(Z)))/n
    ## Variance for each component of random effects
    VarR <-  c(VarSite,
               VarCorr(model)$`spp:site`[1],
               VarCorr(model)$`spp:region`[1],
               VarCorr(model)$`spp`[1])
    names(VarR) <- c("site", "spp:site", "spp:region", "spp")
    ## Denominator for R2GLMM formula works for Poisson distribution only
    deno <- (VarF + sum(VarR) + #pi^2/3
             log(1 + 1/exp(as.numeric(fixef(m0)))))
    ## R2GLMM(m) - marginal R2GLMM 
    r2f <- VarF/deno
    ## R2GLMM(c) - conditional R2GLMM for full model
    r2t <- (VarF + sum(VarR))/deno
                                        # R2 random effects only
    r2rand <- r2t-r2f
    ## R2 Residuals
    r2res <- 1-r2t
    ## Partitioning R2 GLMM for each random effect
    r2rand.part <- VarR/deno
    r2.tab <- data.frame(component=c("conditional", "fixed", "random",
                                     names(VarR)),
                         R2=c(r2t,r2f,r2rand, r2rand.part),
                         type=c("all", "niche", "all.random", "niche", "neutral", "neutral", "idiosyncratic"))
    r2.partition <- aggregate(r2.tab$R2, list(type=r2.tab$type), sum)
    names(r2.partition)[2] <- "R2.partition"
    R2 <- r2.tab$R2
    names(R2) <- r2.tab$component
    return(list(full.table=r2.tab, part.table=r2.partition, R2=R2))
}

## Function for full2 model
## added 'null.model' argument and '...' argument to the fit of the null model (to allow optmization control)
r2.full2 <- function(model, null.model, ...){
    ## Function to calculate the null model, model with all random terms 
    best.null <- function(model) {
        parens <- function(x) paste0("(",x,")")
        onlyBars <- function(form) reformulate(sapply(findbars(form),
                                                      function(x)  parens(deparse(x))),
                                               response=".")
        onlyBars(formula(model))
        best.null <- update(model,onlyBars(formula(model)), ...)
        return(best.null)
    }
    ## Calculates null model
    if(missing(null.model))
        m0 <- best.null(model)
    else
        m0 <- null.model
    ## Variance for fixed effects
    VarF <- var(as.vector(fixef(model) %*% t(model@pp$X)))
    ## Variance for slope/intercept random effects
    ## site    
    X <- model.matrix(model)
    n <- nrow(X)
    Z <- X[, c("(Intercept)","grad")]
    sigma.site <- VarCorr(model)$site
    VarSite <- sum(diag(Z %*% sigma.site %*% t(Z)))/n
    ## Variance for slope/intercept random effects
    ## spp
    sigma.spp <- VarCorr(model)$spp
    VarSpp <- sum(diag(Z %*% sigma.spp %*% t(Z)))/n
    ## Variance for each component of random effects
    VarR <-  c(VarSite, VarSpp,
               VarCorr(model)$`spp:site`[1],
               VarCorr(model)$`spp:region`[1])
    names(VarR) <- c("site", "spp", "spp:site", "spp:region")
                                        # Denominator for R2GLMM formula works for Poisson distribution only
    deno <- (VarF + sum(VarR) + #pi^2/3
             log(1 + 1/exp(as.numeric(fixef(m0)))))
                                        # R2GLMM(m) - marginal R2GLMM 
    r2f <- VarF/deno
                                        # R2GLMM(c) - conditional R2GLMM for full model
    r2t <- (VarF + sum(VarR))/deno
                                        # R2 random effects only
    r2rand <- r2t-r2f
    ## R2 Residuals
    r2res <- 1-r2t
    ## Partitioning R2 GLMM for each random effect
    r2rand.part <- VarR/deno
    r2.tab <- data.frame(component=c("conditional", "fixed", "random",
                                     names(VarR)),
                         R2=c(r2t,r2f,r2rand, r2rand.part),
                         type=c("all", "niche", "all.random", "niche",  "niche", "neutral", "neutral"))
    r2.partition <- aggregate(r2.tab$R2, list(type=r2.tab$type), sum)
    names(r2.partition)[2] <- "R2.partition"
    R2 <- r2.tab$R2
    names(R2) <- r2.tab$component
    return(list(full.table=r2.tab, part.table=r2.partition, R2=R2))
}

## Funtion for neutral model
## added 'null.model' argument and '...' argument to the fit of the null model (to allow optmization control)
r2.neutral <- function(model, null.model, ...){
    ## Function to calculate the null model, model with all random terms 
    best.null <- function(model) {
        parens <- function(x) paste0("(",x,")")
        onlyBars <- function(form) reformulate(sapply(findbars(form),
                                                      function(x)  parens(deparse(x))),
                                               response=".")
        onlyBars(formula(model))
        best.null <- update(model,onlyBars(formula(model)), ...)
        return(best.null)
    }
    ## Calculates null model
    if(missing(null.model))
       m0 <- best.null(model)
    else
        m0 <- null.model
    ## Variance for fixed effects
    VarF <- var(as.vector(fixef(model) %*% t(model@pp$X)))
    ## Variance for each component of random effects
    VarR <-  c(VarCorr(model)$`site`[1],
               VarCorr(model)$`spp:site`[1],
               VarCorr(model)$`spp:region`[1],
               VarCorr(model)$`spp`[1]
               )
    names(VarR) <- c("site", "spp:site", "spp:region", "spp")
    ## Denominator for R2GLMM formula works for Poisson distribution only
    deno <- (VarF + sum(VarR) + #pi^2/3
             log(1 + 1/exp(as.numeric(fixef(m0)))))
    ## R2GLMM(m) - marginal R2GLMM 
    r2f <- VarF/deno
    ## R2GLMM(c) - conditional R2GLMM for full model
    r2t <- (VarF + sum(VarR))/deno
    ## R2 random effects only
    r2rand <- r2t-r2f
    ## R2 Residuals
    r2res <- 1-r2t
    ## Partitioning R2 GLMM for each random effect
    r2rand.part <- VarR/deno
    r2.tab <- data.frame(component=c("conditional", "fixed", "random",
                                     names(VarR)),
                         R2=c(r2t,r2f,r2rand, r2rand.part),
                         type=c("all", "niche", "all.random", "neutral",  "neutral", "neutral", "idiosyncratic"))
    r2.partition <- aggregate(r2.tab$R2, list(type=r2.tab$type), sum)
    names(r2.partition)[2] <- "R2.partition"
    R2 <- r2.tab$R2
    names(R2) <- r2.tab$component
    return(list(full.table=r2.tab, part.table=r2.partition, R2=R2))
}

