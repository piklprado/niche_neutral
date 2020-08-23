#####################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al #####################################################################
#####################################################################################

## Code for applying the niche-neutral GLMM framework
## Data and analysis from the manuscript

###############################
# PART 1: loading packages ####
###############################
# required packages
library(bbmle)
library(lme4)
# library(optimx)
library(xtable)
library(piecewiseSEM)
library(dplyr)
library(ggplot2)
library(sads)
# functions to calculate R2
source("code/functions.R")

############################
# PART 2: loading data #####
############################

## abundant
fern.data.ab <- read.csv("data/data_ab_for_modeling.csv", as.is = TRUE)

## rare
fern.data.rare <- read.csv("data/data_rare_for_modeling.csv", as.is = TRUE)

## all
fern.data <- read.csv("data/data_for_modeling.csv", as.is = TRUE)

#######################
# PART 3: modeling ####
#######################
## Building models corresponding to general hypothesis

### Fitting models only for abundant species
head(fern.data.ab)
head(fern.data.rare)

## 3.1: abundant species

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m.full.ab <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                   ##+ indumentum*grad + indumentum*I(grad^2)
                   +  life_form*grad + life_form*I(grad^2)
                   + (1|spp) + (1|spp:region) + (1|spp:site) + (1 + grad|site),
                   data = fern.data.ab, family = "poisson",
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 5e6)))

m.full2.ab <- glmer(abundance ~ grad + I(grad^2)
                    ##+ (1|spp)
                    + (1|spp:region) + (1|spp:site) + (1 + grad|spp) + (1 + grad|site),
                    data = fern.data.ab, family = "poisson",
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 1e6)))

# ### using estimated parameters to buit a new model
ss <- getME(m.full2.ab, c("theta","fixef"))
m.full2b.ab <- update(m.full2.ab, start = ss)

m.neutral.ab <- glmer(abundance ~ (1|spp:region) + (1|spp:site) + (1|site),
                      data = fern.data.ab, family = "poisson",
                      control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 1e6)))

m.niche.ab <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                    ##+ indumentum*grad + indumentum*I(grad^2)
                    +  life_form*grad + life_form*I(grad^2)
                    + (1|spp) + (1 + grad|site),
                    data = fern.data.ab, family = "poisson",
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 5e6)))

m.env.ab <- glmer(abundance ~ grad + I(grad^2)
                  + (1+grad|site) + (1+grad|spp),
                  data = fern.data.ab, family = "poisson",
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 5e6)))

m.null.ab <- glmer(abundance ~ 1
                   + (1|spp) + (1|region) + (1|site),
                   data = fern.data.ab, family = "poisson",
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 5e6)))

## Model selection
m.list <- list(m.full.ab, m.neutral.ab, m.niche.ab, m.null.ab, m.full2b.ab,
               m.env.ab)#, m.full.thickness, m.full.lifeform, m.full.indument)
mod.names <- c("niche & neutral", "neutral", "niche", "null", "env & neutral",
               "env")#, "thickness&neutral", "lifef&neutral", "indument&neutral")

### AIC
aic.ab <- AICctab(m.list, mnames = mod.names, base = TRUE, weights = TRUE, logLik = TRUE)
class(aic.ab) <- 'data.frame'

write.table(aic.ab, "../results/aic_table_ab.csv",
            row.names = TRUE, col.names = TRUE, sep = ",")

# 3.2 Rare species

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among spp sharing the same ES, local and regional limited dispersal
m.full.rar <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                    ##+ indumentum*grad + indumentum*I(grad^2)
                    +  life_form*grad + life_form*I(grad^2)
                    + (1|spp) + (1|spp:region) + (1|spp:site) + (1 + grad|site),
                    data = fern.data.rare, family = "poisson",
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 5e6)))

## Trying different combinations of traits
m.full2.rar <- glmer(abundance ~ grad + I(grad^2)
                     ##+ (1|spp)
                     + (1|spp:region) + (1|spp:site) + (1+grad|spp) + (1 + grad|site),
                     data = fern.data.rare, family = "poisson",
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 5e6)))

m.neutral.rar <- glmer(abundance ~ (1|spp:region) + (1|spp:site) + (1|site),
                       data = fern.data.rare, family = "poisson",
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 1e6)))

m.niche.rar <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                     ##+ indumentum*grad + indumentum*I(grad^2)
                     +  life_form*grad + life_form*I(grad^2)
                     + (1|spp) + (1 + grad|site),
                     data = fern.data.rare, family = "poisson",
                     control = glmerControl(optimizer = "bobyqa",
                                            optCtrl = list(maxfun = 5e6)))

m.env.rar <- glmer(abundance ~ grad + I(grad^2)
                   + (1 + grad|site) + (1 + grad|spp),
                   data = fern.data.rare, family = "poisson",
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 5e6)))

m.null.rar <- glmer(abundance ~ 1
                    + (1|spp) + (1|region) + (1|site),
                    data = fern.data.rare, family = "poisson",
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 5e6)))

## Model selection
m.list.rar <- list(m.full.rar, m.neutral.rar, m.niche.rar, m.null.rar, m.full2.rar,
               m.env.rar)#, m.full.thickness, m.full.lifeform, m.full.indument)

### AIC
aic.rar <- AICctab(m.list.rar, mnames = mod.names, base = TRUE, weights = TRUE,
                   logLik = TRUE)
class(aic.rar) <- 'data.frame'

write.table(aic.rar, "../results/aic_table_rar.csv",
            row.names = TRUE, col.names = TRUE, sep = ",")


##############################################
## Part 4 : calculating R2 ###################
##############################################

### using functions from tutorial in file functions.R
r2.ab <- r2.full(m.full.ab)
r2.rar <- r2.neutral(m.neutral.rar)

r2.ab
r2.rar

#### Grafico R2 ####
## creating data frame for plot
r2.df <- rbind(r2.ab$full.table, r2.rar$full.table)
r2.df$group <- c(rep("abundant", 7), rep("rare", 6))

r2.df <- r2.df[!r2.df$type %in% c("all", "all.random"),]

#pdf("../figures/barplot.pdf")
ggplot(r2.df, aes(fill = type, y = R2, x = group)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic()
#dev.off()

##############################################
## Part 5: exporting predicted data ##########
##############################################

## Predicted for abundant species by the selected model
pred.ab <- merTools::predictInterval(m.full.ab,  which = "all",
                                     type = "linear.prediction", level = 0.95)
## Predicted for rare species
pred.rare <- merTools::predictInterval(m.neutral.rar,  which = "all",
                                       type = "linear.prediction", level = 0.95)

## Mean expected abundances in function of altitude for all species by the combination of the best models
all.preds <- rbind(cbind(fern.data.ab,combined = pred.ab[pred.ab$effect == "combined", "fit"],
                         fixed = pred.ab[pred.ab$effect == "fixed", "fit"]),
                   cbind(fern.data.rare, combined = pred.rare[pred.rare$effect == "combined","fit"],
                         fixed = pred.rare[pred.rare$effect == "fixed","fit"]))

## Log of mean of observed values, log of sds of predicted by total and random effects + intercepts
sads.meta <- all.preds %>%
    mutate(random = combined - fixed +
               ifelse(ab.rare == "rare", fixef(m.neutral.rar)[1] , fixef(m.full.ab)[1]),
           ab.class=factor(ifelse(ab.rare == "abundant", "Core", "Occasional"))) %>%
    group_by(spp, ab.class) %>%
    summarise(mean.obs = mean(abundance), sd.obs = sd(abundance),
              lwr.obs = ifelse(mean.obs - sd.obs > 0, mean.obs - sd.obs, 1/30),
              upr.obs = mean.obs + sd.obs,
              mean.comb = mean(exp(combined)), sd.comb = sd(exp(combined)),
              mean.random = mean(exp(random)), sd.random = sd(exp(random)),
              lwr.random = ifelse(mean.obs - sd.random > 0, mean.obs - sd.random, 1 / 30),
              upr.random = mean.obs+sd.random) %>%
    ungroup() %>%
    mutate(sp.rank = rank(-mean.obs, ties.method = "first"),
           sp.rank.c = rank(-mean.comb, ties.method = "first"),
           sp.rank.r = rank(-mean.random, ties.method = "first"))


## writing predicted values
write.table(all.preds, 'results/predicted.csv',
            col.names = TRUE, row.names = FALSE, sep = ",")


write.table(sads.meta, 'results/sads_metacommunity.csv',
            col.names = TRUE, row.names = FALSE, sep = ",")
