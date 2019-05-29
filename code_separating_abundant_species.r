######################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al ######################################################################
#####################################################################################

## Code for applying the niche-neutral GLMM framework
## Data and analysis from the manuscript  
### Core for separating abundant and rare spp

################################
# PART 1: loading packages ####
###############################

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
source("functions.R")


#############################
# PART 2: loading data ######
############################

# reading data
fern.data <- read.csv("fern_data_Mortaraetal.csv", as.is=TRUE)
head(fern.data)
# changing column names
names(fern.data)[c(2, 4, 9)] <- c("spp", "region", "grad")
# removing NAs
fern.data <- na.omit(fern.data) # tem 60 NA's de indumento e espessura de folha. É isto mesmo? sao duas especies apenas, mas como todo mundo ocorre em todo lugar ficam 60
#fern.data$site <- scale(rep(1:30, length(unique(fern.data$species)))) ## troquei por um fator
# creating vector with site id
fern.data$site <- factor(rep(1:30, length(unique(fern.data$spp))))
sp <- unique(fern.data$spp)
    
### Creating sp site matrix
sp.site <- cast(fern.data, site + altitude ~ spp, value='abundance', FUN=mean)
sp.site

site <- sp.site$site
altitude <- sp.site$altitude

sp.site <- apply(sp.site[,c(-1, -2)], 2, as.numeric)

dim(sp.site)

sp.alt <- apply(sp.site, 2, function(x) tapply(x, altitude, mean))

## Creating vector of species abundances for the total metacommunity
ab.meta <- colSums(sp.site)
names(ab.meta) <- sp

length(ab.meta)

ab.meta

#### separating 40 most abundant species in the metacommunity
indice.sp <- order(ab.meta, decreasing=TRUE)[1:40]

freq <- rowSums(sp.site>0)

ab.mean <- apply(sp.site, 1, mean)
ab.mean

plot(freq ~ ab.mean)

#### Will indentify species as abundant (first 40 in rank) and rare (other)
ab.rare <- data.frame(spp=unique(fern.data$spp),
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

### Fitting models only for abundant species
head(fern.data.ab)
head(fern.data.rare)

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m.full.ab <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                ##+ indumentum*grad + indumentum*I(grad^2)
                +  life_form*grad + life_form*I(grad^2)
                + (1|spp) + (1|spp:region) + (1|spp:site) + (1+grad|site),
                data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.full2.ab <- glmer(abundance ~ grad + I(grad^2)
                 ##+ (1|spp)
                 + (1|spp:region) + (1|spp:site) + (1+grad|spp) + (1+grad|site),
                  data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

# ### using estimated parameters to buit a new model
ss <- getME(m.full2.ab, c("theta","fixef"))
m.full2b.ab <- update(m.full2.ab, start=ss)

m.neutral.ab <- glmer(abundance ~ (1|spp) + (1|spp:region) + (1|spp:site) + (1|site),
                   data=fern.data.ab, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

m.niche.ab <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                 ##+ indumentum*grad + indumentum*I(grad^2)
                 +  life_form*grad + life_form*I(grad^2)
                 + (1|spp) + (1+grad|site),
                 data=fern.data.ab, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.env.ab <- glmer(abundance ~ grad + I(grad^2)
               + (1+grad|site) + (1+grad|spp),
               data=fern.data.ab, family="poisson",
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.null.ab <- glmer(abundance ~ 1 
                + (1|spp) + (1|region) + (1|site),
                data=fern.data.ab, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Model selection
m.list <- list(m.full.ab, m.neutral.ab, m.niche.ab, m.null.ab, m.full2b.ab,
               m.env.ab)#, m.full.thickness, m.full.lifeform, m.full.indument)
mod.names <- c("niche & neutral", "neutral", "niche", "null", "env & neutral",
               "env")#, "thickness&neutral", "lifef&neutral", "indument&neutral")

### AIC
AICctab(m.list, mnames=mod.names, base=TRUE, weights=TRUE, logLik=TRUE)

#######################
### Rare species now #####
########################

## Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among spp sharing the same ES, local and regional limited dispersal
m.full.rar <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                ##+ indumentum*grad + indumentum*I(grad^2)
                +  life_form*grad + life_form*I(grad^2)
                + (1|spp) + (1|spp:region) + (1|spp:site) + (1+grad|site),
                data=fern.data.rare, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Trying different combinations of traits
m.full2.rar <- glmer(abundance ~ grad + I(grad^2)
                 ##+ (1|spp)
                 + (1|spp:region) + (1|spp:site) + (1+grad|spp) + (1+grad|site),
                 data=fern.data.rare, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.neutral.rar <- glmer(abundance ~ (1|spp) + (1|spp:region) + (1|spp:site) + (1|site),
                   data=fern.data.rare, family="poisson",
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

m.niche.rar <- glmer(abundance ~ thickness*grad + thickness*I(grad^2)
                 ##+ indumentum*grad + indumentum*I(grad^2)
                 +  life_form*grad + life_form*I(grad^2)
                 + (1|spp) + (1+grad|site),
                 data=fern.data.rare, family="poisson",
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.env.rar <- glmer(abundance ~ grad + I(grad^2)
               + (1+grad|site) + (1+grad|spp),
               data=fern.data.rare, family="poisson",
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.null.rar <- glmer(abundance ~ 1 
                + (1|spp) + (1|region) + (1|site),
                data=fern.data.rare, family="poisson",
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

## Model selection
m.list.rar <- list(m.full.rar, m.neutral.rar, m.niche.rar, m.null.rar, m.full2.rar,
               m.env.rar)#, m.full.thickness, m.full.lifeform, m.full.indument)

### AIC
AICctab(m.list.rar, mnames=mod.names, base=TRUE, weights=TRUE, logLik=TRUE)

###############################################
## Part 4 : calculating R2 ####################
##############################################

### using functions from tutorial in file functions.R 
r2.ab <- r2.full(m.full.ab)
r2.rar <- r2.neutral(m.neutral.rar)
r2.rar2 <- r2.full(m.full.rar)

r2.ab
r2.rar


#### Grafico R2 ####
## creating data frame for plot
r2.df <- rbind(r2.ab$full.table, r2.rar$full.table)
r2.df$group <- rep(c("abundant", "rare"), each=7)

r2.df <- r2.df[!r2.df$type%in%c("all", "all.random"),]

ggplot(r2.df, aes(fill=type, y=R2, x=group)) + 
    geom_bar(stat="identity", position="fill") +
    theme_classic()

## Metodo 1: desvios-padrão dos previstos

##################################################################
## Part 5 : generating observed and predicted values ############
##################################################################

### calculando os previstos
pred.values <- predict(m.full.ab, re.form=NULL, type='response') 
# changed to type=response to have abundance back on original scale
pred.values2 <- predict(m.neutral.rar, re.form=NULL, type='response')

# Auxiliar function to calculate mean, sd, lower and upper interval of abundance values
mod.df <- function(values, data) {
mod.table <- aggregate(values, 
                        list(altitude=data[,'altitude'], 
                             thickness=data[,'thickness'],
                             life_form=data[,'life_form']), mean)
mod.table$sd <- aggregate(values, 
                        list(altitude=data[,'altitude'], 
                             thickness=data[,'thickness'],
                             life_form=data[,'life_form']), sd)$x 
mod.table$plwr <- mod.table$x - mod.table$sd
mod.table$pupr <- mod.table$x + mod.table$sd
names(mod.table)[4] <- "Abundance"
return(mod.table)
}

# Predicted values
## table for abundant
pred.table <- mod.df(pred.values, fern.data.ab)
## table for rare 
pred.table2 <- mod.df(pred.values2, fern.data.rare)

# Observed values
## table for abundant
obs <- mod.df(fern.data.ab$abundance, fern.data.ab)
## table for rare 
obs.rar <- mod.df(fern.data.rare$abundance, fern.data.rare)


#############################################
######### Creating figures ##################
############################################

## Defining colours for graphics
cor1 <-rgb(140, 1, 28, maxColorValue=255) #rgb(44, 152, 32, maxColorValue=255) # terrestre
cor3 <- rgb(4, 70, 120, maxColorValue=255) #rgb(239, 144, 33, maxColorValue=255) # hemi
cor2 <- rgb(199, 172, 29, maxColorValue=255) # ep

## Vector for abundant and rare species
cores <- rep(c(cor1, cor2, cor3), each=20)
cores.rar <- c(rep(cor1, 20), rep(cor2, 10), rep(cor3, 20))
## Names for graphics
nomes <- c(ep="epiphyte", hemi="hemiepiphyte", ter="terrestrial", coriacea="coriaceous", membranacea="membranaceous")

#############################################
#### Grafico com a previsao dos modelos ####
#############################################

## Um grafico rapido
#pdf("abudant_gradient.pdf")
obs %>%
    mutate(Altitude=rep(unique(altitude), 6)) %>%
    ggplot(aes(Altitude, Abundance)) +
    ylab("Abundance (log)") + #scale_y_log10() + # making only y axis in log and not abundance values
    facet_grid(thickness ~ life_form, labeller=as_labeller(nomes)) +
    geom_point(colour=cores, size=3) +
    geom_linerange(aes(ymin=plwr, ymax=pupr), colour=cores) +
    geom_ribbon(aes(x=altitude, y=Abundance, ymin=plwr, ymax=pupr), 
                data=pred.table, alpha=0.1) +
    theme_classic(base_size=15) #+
#    theme(strip.background = element_blank(), strip.text = element_blank())
#dev.off()


#########################################
#### Grafico SADS  #####################
#########################################

# creating data frame with predicted mean and sp for each species in the metacommunity
meta.df <- aggregate(fern.data$abundance, list(spp=fern.data$spp), mean)
meta.df$sd.obs <- aggregate(fern.data$abundance, list(fern.data$spp), sd)$x
meta.df$rarity <- ab.rare$ab.rare
head(meta.df)
names(meta.df)[2] <- "mean.obs"

head(meta.df)

# mean and sd for abundant and rare species
ab.df <- aggregate(pred.values,
                   list(spp=fern.data.ab$spp), mean)
ab.df$sd.pred <- aggregate(pred.values, 
                           list(spp=fern.data.ab$spp), sd)$x
rare.df <- aggregate(pred.values2, 
                     list(spp=fern.data.rare$spp), mean)
rare.df$sd.pred <- aggregate(pred.values2, 
                             list(spp=fern.data.rare$spp), sd)$x

ab.rare.df <- bind_rows(ab.df, rare.df)
ab.rare.df.or <- ab.rare.df[order(ab.rare.df$spp),] 
dim(ab.rare.df.or)

meta.df$mean.pred <- ab.rare.df.or$x
meta.df$sd.pred <- ab.rare.df.or$sd

# colors from thesis
# \definecolor{niche}{RGB}{202, 0, 32}
# \definecolor{neutral}{RGB}{4,4,160}
# \definecolor{nineu}{RGB}{245,219,58}
# \definecolor{grey}{RGB}{180,177,172}

ab.cor <- rgb(245,219,58, maxColorValue = 255)
rare.cor <- rgb(4,4,160, maxColorValue = 255)

meta.df$cor <- ifelse(meta.df$rarity=="rare", rare.cor, ab.cor)

# write.table(meta.df, "metacommunity_predicted_abundances.csv", 
#             sep=",", row.names=FALSE, col.names = TRUE)

# sad total
sad <- ggplot(aes(reorder(spp, -mean.obs), mean.obs), data=meta.df) +
  ylab("Abundance (log)") + #+ scale_y_log10() + # making only y axis in log and not abundance values
  xlab("Species") +
  geom_point(colour=meta.df$cor, size=3) +
  geom_linerange(aes(ymin=mean.obs-sd.obs, 
                     ymax=mean.obs+sd.obs), colour=meta.df$cor) +
  #facet_grid(ab.tot$ab.rare) +
   geom_ribbon(aes(spp, 
                   ymin=mean.pred-sd.pred, 
                   ymax=mean.pred+sd.pred), alpha=1) +
  theme_classic(base_size=15) #+
#dev.off()

sad

pdf("grafico_sad.pdf")
sad + theme(axis.text.x = element_blank()) + coord_cartesian(xlim=c(0,153), ylim=c(0,110))
dev.off()

### grafico com as oitavas
oc.com <- octav(meta.df$mean.obs)

## para ps valores observados
oc.abra <- tapply(meta.df$mean.obs, list(meta.df$rarity), octav)
oc.df <- rbind(oc.abra[[1]], oc.abra[[2]])

oc.df$abra <- as.vector(c(rep("abundant", 8), rep("rare", 10)))

# agora para os previstos

### preparando os dados
prev.ab <- data.frame(abundance=pred.values, site=fern.data.ab$site, 
                      altitude=fern.data.ab$altitude, spp=fern.data.ab$spp)

prev.ra <- data.frame(abundance=pred.values2, site=fern.data.rare$site, 
                      altitude=fern.data.rare$altitude, spp=fern.data.rare$spp)

all.prev <- list(prev.ab, prev.ra)

sp.site.prev <- lapply(all.prev, function(x) cast(x, site + altitude ~ spp, 
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


    

