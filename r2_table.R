# Function to calculate R-squared for each random effect based on code from Melina Leite and Nakagawa
# Works only for Poisson GLMMs
r2.table <- function(model){
# Function to calculate the null model, model with all random terms 
best.null <- function(model) {
parens <- function(x) paste0("(",x,")")
onlyBars <- function(form) reformulate(sapply(findbars(form),
                                              function(x)  parens(deparse(x))),
                                              response=".")
onlyBars(formula(model))
best.null <- update(model,onlyBars(formula(model)))
return(best.null)
}
# Calculates null model
m0 <- best.null(model)
# Variance for fixed effects
VarF <- var(as.vector(fixef(model) %*% t(model@pp$X)))
    ## Denominator for R2GLMM formula works for Poisson distribution only
    s2d <- log(1 + 1/exp(as.numeric(fixef(m0)))) ## Este é termo sigma^2_d de Nakagawa et al. 2017 Interface. Mas veja eqs 5.7-5.8 para correcao.
    sigma.tau <- sum(unlist(VarCorr(model))) ## de acordo eq 5.8 Nakagawa et al. 2017
deno <- VarF + sigma.tau +  s2d 
# R2GLMM(m) - marginal R2GLMM 
r2f <- VarF/deno
# R2GLMM(c) - conditional R2GLMM for full model
r2t <- (VarF + sum(unlist(VarCorr(model))))/deno
# R2 random effects only
r2rand <- r2t-r2f
## R2 Residuals
r2res <- 1-r2t
## Partitioning R2 GLMM for each random effect
r2rand.part <- unlist(VarCorr(model))/deno
r2.tab <- t(as.data.frame(c(conditional = r2t,
      fixed = r2f,
      random = r2rand,
      r2rand.part)))
row.names(r2.tab) <- "model"
return(as.data.frame(r2.tab))
}
