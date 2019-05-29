### Funtion for the full model
r2_full <- function(model, ...){
  # Function to calculate the null model, model with all random terms 
  best.null <- function(model, ...) {
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
  X <- model.matrix(model)
  n <- nrow(X)
  Beta <- fixef(model)
  VarF <- var(X %*% Beta)
  ## Variance for slope/intercept random effects
  #Z <- X[, c("(Intercept)","grad")]
  #sigma <- VarCorr(model)$site
  Sigma.list <- VarCorr(model)
  VarR <- sapply(Sigma.list,
             function(Sigma)
             {
               Z <-X[,rownames(Sigma)]
               sum(diag(Z %*% Sigma %*% t(Z)))/n
             })
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
  r2rand.part <- as.vector(VarR)/as.vector(deno)
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
