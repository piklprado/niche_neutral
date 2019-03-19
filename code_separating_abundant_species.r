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
                                               life_form=fern.data.ab$life_form), sd))$x
head(pred.table)
pred.table$plwr <- pred.table$x - pred.table$sd
pred.table$pupr <- pred.table$x + pred.table$sd

### PARA RARAS
# previstos sobre os dados originais (omitindo argumento newdata, usa-se a tabela original de dados)
pred.values2 <- predict(m.neutral.rar, re.form=NULL, type='response') # changed to type=response to have abundance back on original scale

# We calculate mean predicted values and standard error for each altitude 
## Predicted mean values
pred.table2 <- aggregate(pred.values2, list(altitude=fern.data.rare$altitude, thickness=fern.data.rare$thickness,
                                          life_form=fern.data.rare$life_form), mean)
## Predicted stardard deviations
pred.table2$sd <- aggregate(pred.values2, list(altitude=fern.data.rare$altitude, thickness=fern.data.rare$thickness,
                                             life_form=fern.data.rare$life_form), sd))$x
head(pred.table2)
pred.table2$plwr <- pred.table2$x - pred.table2$sd
pred.table2$pupr <- pred.table2$x + pred.table2$sd

### PARA ABUNDANTES
# Second we create a data frame with observed mean values and its standard error
obs <- aggregate(fern.data.ab$abundance, by=list(altitude=fern.data.ab$alt_std, thickness=fern.data.ab$thickness,
                                              life_form=fern.data.ab$life_form), mean)
## Observed standard deviation
obs$std <- aggregate(fern.data.ab$abundance, by=list(altitude=fern.data.ab$alt_std, thickness=fern.data.ab$thickness,
                                              life_form=fern.data.ab$life_form), sd))$x
head(obs)
names(obs) <- c("Altitude", "thickness", "life_form", "Abundance", "std")

### PARA RARAS
# Second we create a data frame with observed mean values and its standard error
obs.rar <- aggregate(fern.data.rare$abundance, by=list(altitude=fern.data.rare$alt_std, thickness=fern.data.rare$thickness,
                                                 life_form=fern.data.rare$life_form), mean)
## Observed standard deviation
obs.rar$std <- aggregate(fern.data.rare$abundance, by=list(altitude=fern.data.rare$alt_std, thickness=fern.data.rare$thickness,
                                                     life_form=fern.data.rare$life_form), sd))$x
head(obs.rar)
names(obs.rar) <- c("Altitude", "thickness", "life_form", "Abundance", "std")

summary(obs)

summary(pred.table)

head(obs)
head(pred.table)

cores <- rep(c(cor1, cor2, cor3), each=20)
cores.rar <- c(rep(cor1, 20), rep(cor2, 10), rep(cor3, 20))
nomes <- c(ep="epiphyte", hemi="hemiepiphyte", ter="terrestrial", coriacea="coriaceous", membranacea="membranaceous")

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

################################################################################
## Metodo 2: com bootMER (muito demorado, desisti)


## Function to run at each boostrap simulation
## library(parallel)
## f1 <- function(., ...)
##     predict(. , re.form = NULL, ...)
## m.full.boot <- bootMer(m.full, f1, nsim=100, parallel="multicore", ncpus=4)
## save.image()
## m.full.boot.fix <- bootMer(m.full, f1, nsim=100, use.u=TRUE, parallel="multicore", ncpus=4)
## save.image()


#############################################
######### Creating figures ##################
############################################

###############################################
######### DAQUI PRA FRENTE AINDA NAO FUNFA ####
###############################################

#########################################
#### GRAFICO SADS  #####################
#########################################

head(com.rank2)

head(atri)

atri.cor <- atri[,c(1, 26, 27)]
head(atri.cor)

atri.cor

atri.cor$comb <- NA
atri.cor$comb2 <- NA


head(atri.cor)
atri.cor$comb[atri.cor$habitoB=="ter" & atri.cor$espessuraB=="membranacea"] <- cor1
atri.cor$comb[atri.cor$habitoB=="ep" & atri.cor$espessuraB=="membranacea"] <- cor3
atri.cor$comb[atri.cor$habitoB=="hemi" & atri.cor$espessuraB=="membranacea"] <- cor2
atri.cor$comb[atri.cor$habitoB=="ter" & atri.cor$espessuraB=="coriacea"] <- cor1
atri.cor$comb[atri.cor$habitoB=="ep" & atri.cor$espessuraB=="coriacea"] <- cor3
atri.cor$comb[atri.cor$habitoB=="hemi" & atri.cor$espessuraB=="coriacea"] <- cor2

atri.cor$comb2[atri.cor$habitoB=="ter" & atri.cor$espessuraB=="membranacea"] <- 1
atri.cor$comb2[atri.cor$habitoB=="ep" & atri.cor$espessuraB=="membranacea"] <- 1
atri.cor$comb2[atri.cor$habitoB=="hemi" & atri.cor$espessuraB=="membranacea"] <- 1
atri.cor$comb2[atri.cor$habitoB=="ter" & atri.cor$espessuraB=="coriacea"] <- 19
atri.cor$comb2[atri.cor$habitoB=="ep" & atri.cor$espessuraB=="coriacea"] <- 19
atri.cor$comb2[atri.cor$habitoB=="hemi" & atri.cor$espessuraB=="coriacea"] <- 19

# funcao para plot das sads com abundancias relativas e atributos
cont.y <- c(1,4,7,10)
cont.x <- 8:10

graf.sad <- function(com=com.rank2, cor=atri.cor$comb, ponto=atri.cor$comb2){
par(mai=c(0.24, 0.6, 0.24, 0.05), oma=c(3, 3, 0.2, 0.1))
layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11,0), 4, 3, byrow=TRUE))
for(i in 1:10){
plot(com.rank2[[i]], log="y", ylim=c(0.0004,0.5), xlim=c(0, 63),
     col=cor[order(com.cota[i,],decreasing=TRUE )][1:riq.cota[i]],
     pch=ponto[order(com.cota[i,],decreasing=TRUE )][1:riq.cota[i]],
     bty="l", cex=1.9, cex.axis=1.5, xlab="", ylab="", las=1,
     yaxt="n", xaxt="n")
mtext(paste(LETTERS[1:10][i], paste(unique(cota)[i], "m", sep=" "), sep=". "), adj=0.05, padj=-0.5, cex=1.2, font=2)
if(i %in% cont.y){
axis(2, las=1, cex.axis=1.5, at=c(0.0005, 0.002, 0.01, 0.05, 0.2), labels=c("0.0005", "0.002", "0.01", "0.05", "0.2"))
}
else{axis(2, at=c(0.0005, 0.002, 0.01, 0.05, 0.2), labels=rep(" ", 5))}
if(i %in% cont.x){
axis(1, las=1, cex.axis=1.5) }
else{axis(1, labels=FALSE)}
}
plot(0,0, axes=FALSE, xlab="", ylab="", col=0)
legend(x=-1.155, y=0.7, c("terrestrial and membranaceous", "terrestrial and coriaceous",
                     "hemiepiphyte and membranaceous", "hemiepiphyte and coriaceous",
                     "epiphyte and membranaceous", "epiphyte and coriaceous"),
                 pch=rep(c(1, 19), 3), col=rep(c(cor1, cor3, cor2), each=2), cex=1.5, pt.cex=1.6, bty="n")
mtext("Species Rank", 1, outer=TRUE, cex=1.3, padj=1)
mtext("Species Relative Abundances (log)", 2, outer=TRUE, cex=1.3, padj=-1)
}




save.image("Mortaraetal.RData")


    

