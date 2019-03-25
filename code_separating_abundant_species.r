######################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al ######################################################################
#####################################################################################

## Code for applying the niche-neutral GLMM framework
## Data and analysis from the manuscript  
### Core for separating abundant and rare species

#############################
# PART 1: loading packages #
############################

# required packages
library(bbmle)
library(lme4)
#library(optimx)
library(xtable)
library(piecewiseSEM)
library(dplyr)
library(ggplot2)
library(reshape)
library(sads)
## functions to calculate R2 
source("r2_full.R")
source("r2_neutral.R")

#############################
# PART 2: loading data ######
############################

fern.data <- read.csv("fern_data_Mortaraetal.csv", header=TRUE)
head(fern.data)
tail(fern.data)
fern.data <- na.omit(fern.data) # tem 60 NA's de indumento e espessura de folha. É isto mesmo? sao duas especies apenas, mas como todo mundo ocorre em todo lugar ficam 60
#fern.data$site <- scale(rep(1:30, length(unique(fern.data$species)))) ## troquei por um fator
fern.data$site <- factor(rep(1:30, length(unique(fern.data$species))))
sp <- unique(fern.data$species)
    
### Creating sp site matrix
sp.site <- cast(fern.data, site + altitude ~ species, value='abundance', FUN=mean)
sp.site

site <- sp.site$site
altitude <- sp.site$altitude

sp.site <- apply(sp.site[,c(-1, -2)], 2, as.numeric)

dim(sp.site)

sp.alt <- apply(sp.site, 2, function(x) tapply(x, altitude, mean))

## Creating vector of species abundances for the total metacommunity
ab.meta <- colSums(sp.site)
names(ab.meta) <- sp

#### separating 40 most abundant species in the metacommunity
indice.sp <- order(ab.meta, decreasing=TRUE)[1:40]

freq <- rowSums(sp.site>0)

ab.mean <- apply(sp.site, 1, mean)
ab.mean

plot(freq ~ ab.mean)


#### Will indentify species as abundant (first 40 in rank) and rare (other)
ab.rare <- data.frame(species=unique(fern.data$species),
                      ab.rare=ifelse(sp%in%sp[indice.sp], "abundant", "rare"))

ab.rare

fern.data.new <- merge(fern.data, ab.rare)
dim(fern.data)
dim(fern.data.new)

head(fern.data.new)

fern.data.ab <- fern.data.new[fern.data.new$ab.rare=="abundant",]
fern.data.rare <- fern.data.new[fern.data.new$ab.rare=="rare",]

#####################################################################
# PART 3: building the model to represent our hypothesis ############
# Step by step building models corresponding to general hypothesis #
####################################################################

#####
### Adjusting models only for abundant species

head(fern.data.ab)
head(fern.data.rare)

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m.full.ab <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
                ##+ indumentum*alt_std + indumentum*I(alt_std^2)
                +  life_form*alt_std + life_form*I(alt_std^2)
                + (1|species) + (1|species:mountain) + (1|species:site) + (1+alt_std|site),
                data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.full2.ab <- glmer(abundance ~ alt_std + I(alt_std^2)
                 ##+ (1|species)
                 + (1|species:mountain) + (1|species:site) + (1+alt_std|species) + (1+alt_std|site),
                  data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

# ### using estimated parameters to buit a new model
# ss <- getME(m.full2.ab, c("theta","fixef"))
# m.full2b.ab <- update(m.full2.ab, start=ss)

m.neutral.ab <- glmer(abundance ~ (1|species) + (1|species:mountain) + (1|species:site) ,#+ (1|site),
                   data=fern.data.ab, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

m.niche.ab <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
                 ##+ indumentum*alt_std + indumentum*I(alt_std^2)
                 +  life_form*alt_std + life_form*I(alt_std^2)
                 + (1|species) + (1+alt_std|site),
                 data=fern.data.ab, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.env.ab <- glmer(abundance ~ alt_std + I(alt_std^2)
               + (1+alt_std|site) + (1+alt_std|species),
               data=fern.data.ab, family="poisson",
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.null.ab <- glmer(abundance ~ 1 
                + (1|species) + (1|mountain) + (1|site),
                data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Model selection
m.list <- list(m.full.ab, m.neutral.ab, m.niche.ab, m.null.ab, m.full2.ab,
               m.env.ab)#, m.full.thickness, m.full.lifeform, m.full.indument)
mod.names <- c("niche & neutral", "neutral", "niche", "null", "env & neutral",
               "env")#, "thickness&neutral", "lifef&neutral", "indument&neutral")

### AIC
AICctab(m.list, mnames=mod.names, base=TRUE, weights=TRUE, logLik=TRUE)
## BIC
BICtab(m.list, mnames=mod.names,base=TRUE, weights=TRUE, logLik=TRUE)

head(fern.data.ab)

#######################
### Rare species now
########################

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m.full.rar <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
                ##+ indumentum*alt_std + indumentum*I(alt_std^2)
                +  life_form*alt_std + life_form*I(alt_std^2)
                + (1|species) + (1|species:mountain) + (1|species:site) + (1+alt_std|site),
                data=fern.data.rare, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Trying different combinations of traits
m.full2.rar <- glmer(abundance ~ alt_std + I(alt_std^2)
                 ##+ (1|species)
                 + (1|species:mountain) + (1|species:site) + (1+alt_std|species) + (1+alt_std|site),
                 data=fern.data.rare, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.neutral.rar <- glmer(abundance ~ (1|species) + (1|species:mountain) + (1|species:site) ,#+ (1|site),
                   data=fern.data.rare, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

m.niche.rar <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
                 ##+ indumentum*alt_std + indumentum*I(alt_std^2)
                 +  life_form*alt_std + life_form*I(alt_std^2)
                 + (1|species) + (1+alt_std|site),
                 data=fern.data.rare, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.env.rar <- glmer(abundance ~ alt_std + I(alt_std^2)
               + (1+alt_std|site) + (1+alt_std|species),
               data=fern.data.rare, family="poisson",
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.null.rar <- glmer(abundance ~ 1 
                + (1|species) + (1|mountain) + (1|site),
                data=fern.data.rare, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Model selection
m.list.rar <- list(m.full.rar, m.neutral.rar, m.niche.rar, m.null.rar, m.full2.rar,
               m.env.rar)#, m.full.thickness, m.full.lifeform, m.full.indument)

### AIC
AICctab(m.list.rar, mnames=mod.names, base=TRUE, weights=TRUE, logLik=TRUE)
## BIC
BICtab(m.list.rar, mnames=mod.names,base=TRUE, weights=TRUE, logLik=TRUE)


###############################################
## Part 4 : calculating R2 ####################
##############################################

### need to fix functions to work here! 
# r2.ab <- r2.full(m.full.ab)
# r2.rar <- r2.neutral(m.neutral.rar)

## Metodo 1: desvios-padrão dos previstos

##################################################################
## Part 5 : generating observed and predicted values ############
##################################################################

### PARA ABUNDANTES
  # previstos sobre os dados originais (omitindo argumento newdata, usa-se a tabela original de dados)
  pred.values <- predict(m.full.ab, re.form=NULL, type='response') # changed to type=response to have abundance back on original scale
  # We calculate mean predicted values and standard error for each altitude 
  ## Predicted mean values
  pred.table <- aggregate(pred.values, list(altitude=fern.data.ab$altitude, thickness=fern.data.ab$thickness,
                                            life_form=fern.data.ab$life_form), mean)
  ## Predicted stardard deviations
  pred.table$sd <- aggregate(pred.values, list(altitude=fern.data.ab$altitude, thickness=fern.data.ab$thickness,
                                               life_form=fern.data.ab$life_form), std)$x
head(pred.table)
pred.table$plwr <- pred.table$x - pred.table$sd
pred.table$pupr <- pred.table$x + pred.table$sd

### PARA RARAS
# previstos sobre os dados originais (omitindo argumento newdata, usa-se a tabela original de dados)
pred.values2 <- predict(m.neutral.rar, re.form=NULL, type='response') # changed to type=response to have abundance back on original scale

## function for standard error
std <- function(x) sd(x)/sqrt(length(x))

# We calculate mean predicted values and standard error for each altitude 
## Predicted mean values
pred.table2 <- aggregate(pred.values2, list(altitude=fern.data.rare$altitude, thickness=fern.data.rare$thickness,
                                          life_form=fern.data.rare$life_form), mean)
## Predicted stardard deviations
pred.table2$sd <- aggregate(pred.values2, list(altitude=fern.data.rare$altitude, thickness=fern.data.rare$thickness,
                                             life_form=fern.data.rare$life_form), std)$x
head(pred.table2)
pred.table2$plwr <- pred.table2$x - pred.table2$sd
pred.table2$pupr <- pred.table2$x + pred.table2$sd

### PARA ABUNDANTES
# Second we create a data frame with observed mean values and its standard error
obs <- aggregate(fern.data.ab$abundance, by=list(altitude=fern.data.ab$alt_std, thickness=fern.data.ab$thickness,
                                              life_form=fern.data.ab$life_form), mean)
## Observed standard deviation
obs$std <- aggregate(fern.data.ab$abundance, by=list(altitude=fern.data.ab$alt_std, thickness=fern.data.ab$thickness,
                                              life_form=fern.data.ab$life_form), std)$x
head(obs)
names(obs) <- c("Altitude", "thickness", "life_form", "Abundance", "std")

### PARA RARAS
# Second we create a data frame with observed mean values and its standard error
obs.rar <- aggregate(fern.data.rare$abundance, by=list(altitude=fern.data.rare$alt_std, thickness=fern.data.rare$thickness,
                                                 life_form=fern.data.rare$life_form), mean)
## Observed standard deviation
obs.rar$std <- aggregate(fern.data.rare$abundance, by=list(altitude=fern.data.rare$alt_std, thickness=fern.data.rare$thickness,
                                                     life_form=fern.data.rare$life_form), std)$x
head(obs.rar)
names(obs.rar) <- c("Altitude", "thickness", "life_form", "Abundance", "std")

summary(obs)

summary(pred.table)

head(obs)
head(pred.table)

#############################################
######### Creating figures ##################
############################################

### cores para os graficos
cor1 <-rgb(140, 1, 28, maxColorValue=255) #rgb(44, 152, 32, maxColorValue=255) # terrestre
cor3 <- rgb(4, 70, 120, maxColorValue=255) #rgb(239, 144, 33, maxColorValue=255) # hemi
cor2 <- rgb(199, 172, 29, maxColorValue=255) # ep


cores <- rep(c(cor1, cor2, cor3), each=20)
cores.rar <- c(rep(cor1, 20), rep(cor2, 10), rep(cor3, 20))
nomes <- c(ep="epiphyte", hemi="hemiepiphyte", ter="terrestrial", coriacea="coriaceous", membranacea="membranaceous")

#############################################
#### Grafico com a previsao dos modelos ####
#############################################

## Um grafico rapido
#pdf("abudant_gradient.pdf")
obs %>%
    mutate(lAb=Abundance, lstd=std, lwr=lAb-lstd, upr=lAb+lstd, 
           Altitude=rep(unique(altitude), 6)) %>%
    ggplot(aes(Altitude, lAb)) +
    ylab("Abundance (log)") + scale_y_log10() + # making only y axis in log and not abundance values
    geom_point(colour=cores, size=3) +
    geom_linerange(aes(ymin=lwr, ymax=upr), colour=cores) +
    facet_grid(thickness ~ life_form, labeller=as_labeller(nomes)) +
    geom_ribbon(aes(x=altitude, y=x, ymin=plwr, ymax=pupr), data=pred.table, alpha=0.1) +
    theme_classic(base_size=15) #+
#    theme(strip.background = element_blank(), strip.text = element_blank())
#dev.off()

#pdf("rare_gradient.pdf")
obs.rar %>%
  mutate(lAb=Abundance, lstd=std, lwr=lAb-lstd, upr=lAb+lstd, 
         Altitude=rep(unique(altitude), 5)) %>%
  ggplot(aes(Altitude, lAb)) +
  ylab("Abundance (log)") + scale_y_log10() + # making only y axis in log and not abundance values
  geom_point(colour=cores.rar, size=3) +
  geom_linerange(aes(ymin=lwr, ymax=upr), colour=cores.rar) +
  facet_grid(thickness ~ life_form, labeller=as_labeller(nomes)) +
  geom_ribbon(aes(x=altitude, y=x, ymin=plwr, ymax=pupr), data=pred.table2, alpha=0.1) +
  theme_classic(base_size=15) #+
#dev.off()

###############################################
######### DAQUI PRA FRENTE AINDA NAO FUNFA ####
###############################################

#########################################
#### Grafico SADS  #####################
#########################################

### abundancia da metacomunidade
ab.tot <- data.frame(species=names(ab.meta), ab=ab.meta)
ab.tot <- merge(ab.tot, ab.rare, by="species")

head(ab.tot)

#ab.med <- lapply(all.data, function(x) x[x$x>0,])

# sad total
sad <- ggplot(aes(reorder(species, -ab), ab), data=ab.tot) +
  ylab("Abundance (log)") + scale_y_log10() + # making only y axis in log and not abundance values
  xlab("Species") +
  geom_point(size=3) +
  #geom_linerange(aes(ymin=lwr, ymax=upr), colour=cores.rar) +
  #facet_grid(ab.tot$ab.rare) +
  #geom_ribbon(aes(x=altitude, y=x, ymin=plwr, ymax=pupr), data=pred.table2, alpha=0.1) +
  theme_classic(base_size=15) #+
#dev.off()

sad

### grafico com as oitavas
oc.com <- octav(ab.tot$ab)

## para ps valores observados
oc.abra <- tapply(ab.tot$ab, list(ab.tot$ab.rare), octav)
oc.df <- rbind(oc.abra[[1]], oc.abra[[2]])

oc.df$abra <- as.vector(c(rep("abundant", 13), rep("rare", 9)))

# agora para os previstos

### preparando os dados
prev.ab <- data.frame(abundance=pred.values, site=fern.data.ab$site, 
                      altitude=fern.data.ab$altitude, species=fern.data.ab$species)

prev.ra <- data.frame(abundance=pred.values2, site=fern.data.rare$site, 
                      altitude=fern.data.rare$altitude, species=fern.data.rare$species)

all.prev <- list(prev.ab, prev.ra)

sp.site.prev <- lapply(all.prev, function(x) cast(x, site + altitude ~ species, 
                                             value='abundance', FUN=mean)) 

### calculando vetor de abundancias previstos
ab.prev <- lapply(sp.site.prev, colSums)
ab.prev

### calculando as oitavas previstas
oc.prev <- lapply(ab.prev, octav)
oc.prev

oc.df.prev <- rbind(oc.prev[[1]], oc.prev[[2]])
oc.df.prev$abra <- c(rep("abundant", 13), rep("rare", 10))

oc.df.prev <- oc.df.prev[,-14]

oc <- ggplot(aes(octave, Freq), data=oc.df) +
  geom_bar(stat="identity") +
  facet_grid(abra ~.) +
  geom_point(aes(octave, Freq), data=oc.df.prev) +
  theme_classic(base_size=15) #+

pdf("octaves.pdf")
oc
dev.off()

#save.image("Mortaraetal.RData")


    

