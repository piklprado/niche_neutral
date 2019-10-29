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
fern.data <- read.csv("../data/fern_data_Mortaraetal.csv", as.is=TRUE)
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
sp.site <- cast(fern.data, site + altitude ~ spp,
                value = 'abundance', FUN = mean)
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

## Heat map of the species x sites matrix
image(x=1:nrow(sp.site), y=1:ncol(sp.site), log(sp.site[,order(-ab.meta)]))
    
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

#write.table(fern.data.new, "../results/fern_data_output.csv",
#            col.names=FALSE, row.names=TRUE, sep=',')

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

m.neutral.ab <- glmer(abundance ~
                                        #    (1|spp)
                          (1|spp:region) + (1|spp:site) + (1|site),
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
### Rare species #####
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

m.neutral.rar <- glmer(abundance ~
                           ##    (1|spp) +
                           (1|spp:region) + (1|spp:site) + (1|site),
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

r2.ab
r2.rar

#### Grafico R2 ####
## creating data frame for plot
r2.df <- rbind(r2.ab$full.table, r2.rar$full.table)
r2.df$group <- rep(c("abundant", "rare"), each=7)

r2.df <- r2.df[!r2.df$type%in%c("all", "all.random"),]

#pdf("../figures/barplot.pdf")
ggplot(r2.df, aes(fill=type, y=R2, x=group)) + 
    geom_bar(stat="identity", position="fill") +
    theme_classic()
#dev.off()

## Metodo 1: desvios-padrão dos previstos

##################################################################
## Part 5 : generating observed and predicted values ############
##################################################################

### calculando os previstos
pred.values <- predict(m.full.ab, re.form=NULL, type='response')
pred.values.log <- predict(m.full.ab, re.form=NULL) 
# changed to type=response to have abundance back on original scale
pred.values2 <- predict(m.neutral.rar, re.form=NULL, type='response')

# Auxiliar function to calculate mean, sd, lower and upper interval of abundance values
mod.df <- function(values, data, log.val=FALSE) {
    if(log.val)
        values <- ifelse(values>0 , log(values), NA)
    mod.table <- aggregate(values, 
                        list(altitude=data[,'altitude'], 
                             thickness=data[,'thickness'],
                             life_form=data[,'life_form']), mean, na.rm=TRUE)
    mod.table$sd <- aggregate(values, 
                              list(altitude=data[,'altitude'], 
                                   thickness=data[,'thickness'],
                                   life_form=data[,'life_form']), sd, na.rm=TRUE)$x 
    mod.table$plwr <- mod.table$x - mod.table$sd
    mod.table$pupr <- mod.table$x + mod.table$sd
    names(mod.table)[4] <- "Abundance"
    return(mod.table)
}

# Predicted values
## table for abundant
pred.table <- mod.df(pred.values, fern.data.ab)
pred.table.log <- mod.df(pred.values.log, fern.data.ab)
## table for rare 
pred.table2 <- mod.df(pred.values2, fern.data.rare)

# Observed values
## table for abundant
obs <- mod.df(fern.data.ab$abundance, fern.data.ab)
obs.log <- mod.df(fern.data.ab$abundance+0.1, fern.data.ab, log=TRUE)
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

## Acho que é melhor usar o predictInterval do que predict para calcular previstos e intervalos
## Mil simulacçoes dos previstos para cada observacao, efeitos fixos e aleatorios
pred.ab.log <- merTools::predictInterval(m.full.ab, which="full", type="linear.prediction", level=0.95)
pred.grad.ab <- cbind(fern.data.ab, pred.ab.log) %>%
    group_by(altitude,thickness,life_form) %>%
    summarise(Abundance=mean(exp(fit)), plwr=quantile(exp(fit),0.025), pupr=quantile(exp(fit),0.975))

## Um grafico rapido
#pdf("../figures/abundant_gradient.pdf")
obs %>%
    ##    mutate(Altitude=rep(unique(altitude), 6)) %>%
    ggplot(aes(altitude, Abundance)) +
    ylab("Abundance") +
    #scale_y_log10() + # making only y axis in log and not abundance values
    facet_wrap(thickness ~ life_form, labeller=as_labeller(nomes), scale="free") +
    geom_point(colour=cores, size=3) +
    geom_linerange(aes(ymin=plwr, ymax=pupr), colour=cores) +
    geom_ribbon(aes(x=altitude, y=Abundance, ymin=plwr, ymax=pupr), 
                data=pred.grad.ab, alpha=0.1)+
    theme_classic(base_size=15) #+
#    theme(strip.background = element_blank(), strip.text = element_blank())
#dev.off()

## Uma alternativa melhor: os pontos jittered, linha do previsto e barras do IC 95% do previsto
fern.data.ab %>%
    ggplot(aes(factor(altitude), abundance+0.1)) +
    ylab("Abundance") +
    scale_y_log10() + # making only y axis in log and not abundance values
    facet_wrap(thickness ~ life_form, labeller=as_labeller(nomes), scale="free") +
    geom_point(position="jitter", col="blue", size=0.75) +
#    geom_boxplot(color="gray", alpha=0.5) + ## boxplots não estão ajudando muito
    geom_line(aes(x=factor(altitude), y=Abundance, group=1), data=pred.grad.ab) +
    geom_linerange(aes(x=factor(altitude), y=Abundance, ymin=plwr, ymax=pupr), 
                data=pred.grad.ab, alpha=0.5) +
    theme_classic(base_size=15)


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
    geom_ribbon(aes(x=reorder(spp, -mean.obs),
                    y=mean.obs, 
                    ymin=mean.pred-sd.pred, 
                    ymax=mean.pred+sd.pred), alpha=0.5) +
  theme_classic(base_size=15) #+
#dev.off()

sad

#pdf("../figures/grafico_sad.pdf")
sad + theme(axis.text.x = element_blank()) + coord_cartesian(xlim=c(0,153), ylim=c(0,110))
#dev.off()

### grafico com as oitavas

### preparando os dados observados
obs.ab <- data.frame(abundance=fern.data.ab$abundance, site=fern.data.ab$site, 
                      altitude=fern.data.ab$altitude, spp=fern.data.ab$spp)

obs.ra <- data.frame(abundance=fern.data.rare$abundance, site=fern.data.rare$site, 
                      altitude=fern.data.rare$altitude, spp=fern.data.rare$spp)

all.obs <- list(obs.ab, obs.ra)

head(all.obs[[1]])

sp.site.obs <- lapply(all.obs, function(x) cast(x, site + altitude ~ spp, 
                                             value='abundance', FUN=mean)) 

### calculando vetor de abundancias previstos
ab.obs <- lapply(sp.site.obs, colSums)
ab.obs

### calculando as oitavas previstas
oc.obs <- lapply(ab.obs, octav)
oc.obs

oc.df.obs <- rbind(oc.obs[[1]], oc.obs[[2]])
oc.df.obs$abra <- c(rep("abundant", 13), rep("rare", 9))

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

oc.df.prev <- oc.df.prev[oc.df.prev$octave>0,]


oc <- ggplot(aes(octave, Freq), data=oc.df.obs) +
  geom_bar(stat="identity") +
  facet_grid(abra ~.) +
  geom_point(aes(octave, Freq), data=oc.df.prev) +
  theme_classic(base_size=15) #+

oc

pdf("../figures/octaves.pdf")
oc
dev.off()

#save.image("Mortaraetal.RData")

##### NEW Code for sads ####

## Predicted for abundant species by the selected model
pred.ab <- merTools::predictInterval(m.full.ab,  which="all", type="linear.prediction", level=0.95)
## Predicted for rare species
pred.rare <- merTools::predictInterval(m.neutral.rar,  which="all", type="linear.prediction", level=0.95)

## Mean expected abundances in function of altitude for all species by the combination of the best models
all.preds <- rbind( cbind(fern.data.ab,combined=pred.ab[pred.ab$effect=="combined","fit"],
                          fixed=pred.ab[pred.ab$effect=="fixed","fit"]),
                   cbind(fern.data.rare, combined=pred.rare[pred.rare$effect=="combined","fit"],
                         fixed =pred.rare[pred.rare$effect=="fixed","fit"]))

## writing predicted values
#write.table(all.preds, '../results/predicted.csv',
#            col.names=TRUE, row.names=FALSE, sep=",")

## Log of mean of observed values, log of sds of predicted by total and random effects + intercepts
sads.meta <- all.preds %>%
    mutate(random = combined - fixed +
               ifelse(ab.rare=="rare", fixef(m.neutral.rar)[1] , fixef(m.full.ab)[1]),
           ab.class=factor(ifelse(ab.rare=="abundant", "Core", "Occasional"))) %>%
    group_by(spp, ab.class) %>%
    summarise(mean.obs=mean(abundance), sd.obs=sd(abundance),
              lwr.obs=ifelse(mean.obs-sd.obs>0,mean.obs-sd.obs,1/30), upr.obs=mean.obs+sd.obs,
              mean.comb = mean(exp(combined)), sd.comb=sd(exp(combined)),
              mean.random = mean(exp(random)), sd.random=sd(exp(random)),
              lwr.random=ifelse(mean.obs-sd.random>0,mean.obs-sd.random,1/30), upr.random=mean.obs+sd.random
              )%>%
    ungroup() %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"),
           sp.rank.c=rank(-mean.comb, ties.method="first"),
           sp.rank.r=rank(-mean.random, ties.method="first"))
### The figure ###
## 1st panel of the figure: RAD with total and random part of standard deviations
fig.meta1 <- sads.meta %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_linerange(aes(x=sp.rank, ymin=lwr.random, ymax=upr.random, color=ab.class), size=0.5) +
    geom_point(aes(color=ab.class)) +
    geom_ribbon(aes(ymin=lwr.obs, ymax=upr.obs), alpha=0.2) +
    labs(x = "Abundance rank", y = "Mean abundance", color = "") +
    scale_y_log10()
## 2nd panel: Abundances of selected species over the Gradient
## Most abundant of core species
fig.meta2 <- fern.data[fern.data$spp==sads.meta$spp[sads.meta$sp.rank==1],]%>%
    ggplot(aes(altitude, abundance)) +
    stat_smooth(se=FALSE, aes(color=region), alpha=0.5, span=0.6) +
    geom_jitter(aes(color=region), size=2.5, alpha=0.5) +
    labs(title = expression(bold("Core:")~italic("Polybotrya cylindrica")),
                            x = "", y = "Abundance", color = "Region")
## 3rd panel: Most abundant of occasional species
fig.meta3 <- fig.meta2 %+% fern.data[fern.data$spp==sads.meta$spp[sads.meta$sp.rank==41],] +
    labs(title = expression(bold("Occasional:")~italic("Trichomanes polypodioides")),
                            x = "", y = "", color = "Region")
## Arranging all panels in a sigle figure, publication quality (see functions.R)
## List of panels
fig.meta.list <- list(
    fig.meta1 + scale_colour_Publication()+ theme_Publication(),
    fig.meta2 + scale_colour_Publication()+ theme_Publication() + theme(legend.position = "none"),
    fig.meta3 + scale_colour_Publication()+ theme_Publication()
)
## Arrange with grid.arrange
cairo_pdf("../figures/rad_mettacomunity.pdf", width = 9, height = 7.5)
grid.arrange(
    grobs=fig.meta.list,
    bottom=textGrob("Altitude (m)", gp=gpar(fontface = "bold", cex=1.2), vjust=-1.25),
    layout_matrix=matrix(c(1,1,2,3), ncol=2, byrow=TRUE)
)
dev.off()



################################################################################
## Old stuff: Some other tries for the plot of SADs
## Some may be usefull for exploration, maybe keep this part in separate code
################################################################################
head(fern.data.ab)

## Predicted for a single region
spp <- unique(fern.data.ab$spp)
grad <- unique(fern.data.ab$grad)
site <- unique(fern.data.ab$site)
region <- unique(fern.data.ab$region)

newdata <- data.frame(spp=rep(spp, each=length(site)),
                      region=rep(region[1], each=length(spp)*length(site)),
                      site=rep(site, each=length(spp)),
                      grad=rep(grad, each=length(spp)))

head(newdata)

mf2.pi.all <- merTools::predictInterval(m.full.ab,  newdata, which="all", type="linear.prediction", level=0.95)

## Predicted abundances, for neutral models for both specie groups
all.preds.n <- rbind( cbind(fern.data.ab,combined=pred.ab.neu[pred.ab.neu$effect=="combined","fit"],
                          fixed=pred.ab.neu[pred.ab.neu$effect=="fixed","fit"]),
                   cbind(fern.data.rare, combined=pred.rare[pred.rare$effect=="combined","fit"],
                         fixed =pred.rare[pred.rare$effect=="fixed","fit"]))
## Prediction by altitude: whole model and only the random part + intercept
all.preds %>%
    mutate(random = combined - fixed + ifelse(ab.rare=="rare", fixef(m.neutral.rar)[1] , fixef(m.full.ab)[1])) %>%
    group_by(spp,altitude, ab.rare) %>%
    summarise(mean.obs=mean(abundance), mean.comb = mean(exp(combined)), mean.random = mean(exp(random))) %>%
    filter(mean.obs>0) %>%
    group_by(altitude) %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"),
           sp.rank.c=rank(-mean.comb, ties.method="first"),
           sp.rank.r=rank(-mean.random, ties.method="first")) %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_point(aes(color=ab.rare)) +
    geom_line(aes(sp.rank.r,mean.random), color="blue") +
    geom_line(aes(sp.rank.c,mean.comb)) +
    scale_y_log10()+
    facet_wrap(~altitude)


## Mean expected abundances in function of altitude for all species by the combination of neutral models
all.preds.n %>%
    group_by(spp,altitude, ab.rare) %>%
    summarise(mean.obs=mean(abundance), mean.pred = mean(exp(combined))) %>%
    filter(mean.obs>0) %>%
    group_by(altitude) %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"), sp.rank.p=rank(-mean.pred, ties.method="first")) %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_point(aes(color=ab.rare)) +
    geom_line(aes(sp.rank,mean.pred)) +
    scale_y_log10()+
    facet_wrap(~altitude)

## Comparing observed x predicted abundances by altitude for the neutral model
all.preds.n %>%
    group_by(spp,altitude, ab.rare) %>%
    summarise(mean.obs=mean(abundance), mean.pred = mean(exp(combined))) %>%
    filter(mean.obs>0) %>%
    group_by(altitude) %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"), sp.rank.p=rank(-mean.pred, ties.method="first")) %>%
    ggplot(aes(mean.obs, mean.pred)) +
    geom_point(aes(color=ab.rare)) +
    geom_abline(intercept=0,slope=1) +
    scale_y_log10()+
    scale_x_log10()+
    facet_wrap(~altitude)

## Average sads at each altitude
fern.data.new %>%
    group_by(spp,altitude, ab.rare) %>%
    summarise(mean.obs=mean(abundance), ab.sd=sd(abundance)) %>%
    filter(mean.obs>0) %>%
    group_by(altitude) %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first")) %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_point(aes(color=ab.rare))+
    scale_y_log10()+
    facet_wrap(~altitude)

## Sads of each region by altitude
fern.data.new %>%
    filter(abundance>0) %>%
    group_by(region,altitude) %>%
    mutate(sp.rank=rank(-abundance, ties.method="first")) %>%
    ggplot(aes(sp.rank, abundance, group=region)) +
    geom_point(aes(color=region))+
    geom_line(aes(color=region))+
    scale_y_log10()+
    facet_wrap(~altitude)
## With predicted by neutral model
all.preds.n %>%
    filter(exp(combined)>=1) %>%
    group_by(region,altitude) %>%
    mutate(sp.rank=rank(-abundance, ties.method="first"), sp.rank.p=rank(-combined, ties.method="first")) %>%
    ggplot(aes(sp.rank, abundance, group=region)) +
    geom_point(aes(color=region))+
    #geom_line(aes(color=region))+
    geom_line(aes(sp.rank.p, exp(combined), color=region)) +
    scale_y_log10()+
    facet_wrap(~altitude)
    

## Metacommunity sad
fern.data.new %>%
    group_by(spp, ab.rare) %>%
    summarise(mean.obs=mean(abundance), ab.sd=sd(abundance)) %>%
    filter(mean.obs>0) %>%
    ungroup() %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first")) %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_point(aes(color=ab.rare))+
    scale_y_log10()
    
## Metacommunity with predicted by neutral and neutral(rare) + full (abundante) model
pred.meta.n <- all.preds.n %>%
    mutate(random = combined - fixed + ifelse(ab.rare=="rare", fixef(m.neutral.rar)[1] , fixef(m.full.ab)[1])) %>%
    group_by(spp, ab.rare) %>%
    summarise(mean.obs=sum(abundance), mean.comb = sum(exp(combined)), mean.random = sum(exp(random))) %>%
    ungroup() %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"),
           sp.rank.c=rank(-mean.comb, ties.method="first"),
           sp.rank.r=rank(-mean.random, ties.method="first"))
## Rank-abundance plot for observed x predicted by neutral (rare) + full (abundance) models
all.preds %>%
    mutate(random = combined - fixed + ifelse(ab.rare=="rare", fixef(m.neutral.rar)[1] , fixef(m.full.ab)[1])) %>%
    group_by(spp, ab.rare) %>%
    summarise(mean.obs=sum(abundance), mean.comb = sum(exp(combined)), mean.random = sum(exp(random))) %>%
    ungroup() %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"),
           sp.rank.c=rank(-mean.comb, ties.method="first"),
           sp.rank.r=rank(-mean.random, ties.method="first")) %>%
    ggplot(aes(sp.rank, mean.obs)) +
    geom_point(aes(color=ab.rare)) +
    ##geom_line(aes(sp.rank.r,mean.random), color="blue") +
    geom_line(aes(sp.rank.c,mean.comb), data=pred.meta.n, col="blue") +
    geom_line(aes(sp.rank.c,mean.comb)) +   
    scale_y_log10()

## Comparing observed x predicted abundances in the metacommunity for the neutral model
all.preds.n %>%
    group_by(spp,ab.rare) %>%
    summarise(mean.obs=mean(abundance), mean.pred = mean(exp(combined))) %>%
    filter(mean.obs>0) %>%
    ungroup() %>%
    mutate(sp.rank=rank(-mean.obs, ties.method="first"), sp.rank.p=rank(-mean.pred, ties.method="first")) %>%
    ggplot(aes(mean.obs, mean.pred)) +
    geom_point(aes(color=ab.rare)) +
    geom_abline(intercept=0,slope=1) +
    scale_y_log10()+
    scale_x_log10()

## Standard deviations by random and combined effects
all.preds %>%
    mutate(random = combined - fixed ) %>%
    group_by(spp, ab.rare) %>%
    summarise(ab.obs=mean(abundance), sd.obs=sd(abundance),
              ab.comb=mean(exp(combined)),  ab.random=mean(exp(random)),
              sd.comb = sd(exp(combined)), sd.random = sd(exp(random))) %>%
    ungroup() %>%
    mutate(sp.rank=rank(-ab.obs, ties.method="first"),
           ab.comb = exp(mu.comb+sig.comb^2/2), ab.random= exp(mu.random+sig.random^2/2))%>%
    ggplot(aes(sp.rank, ab.obs)) +
    geom_point(aes(color=ab.rare)) +
    geom_linerange(aes(ymin=exp(ab.comb-sd.comb), ymax=exp(ab.comb+sd.comb), color=ab.rare)) +
    geom_linerange(aes(ymin=exp(ab.comb-sd.random), ymax=ab.comb+sd.random)) +
    scale_y_log10()

## Mean rads for both models
pred.ab2 <- merTools::predictInterval(m.full.ab,  which="full", type="linear.prediction", level=0.95, returnSims=TRUE)
pred.ab.neu2 <- merTools::predictInterval(m.neutral.ab,  which="full", type="linear.prediction", level=0.95, returnSims=TRUE)
pred.rar2 <- merTools::predictInterval(m.neutral.rar,  which="full", type="linear.prediction", level=0.95, returnSims=TRUE)

tmp1 <- rbind(cbind(fern.data.ab, attr(pred.ab2, "sim.results")),
      cbind(fern.data.rare, attr(pred.rar2, "sim.results")))
tmp2 <- aggregate(exp(tmp1[,12:1011]), by=list(spp=tmp1$spp), sum)
tmp3 <- apply(tmp2[,-1], 2, sort, decreasing=TRUE)
tmp4 <- aggregate(all.preds$abundance, by=list(spp=all.preds$spp), sum)
plot(rad(tmp4$x))
lines(rad(apply(tmp3, 1, mean)))

plot(rad(tmp3[,1]), ylim=range(tmp3), type="n")
for(i in sample(ncol(tmp3),100))
    lines(rad(tmp3[,i]), col="grey")
lines(rad(tmp4$x))

## Apenas abundantes
tmp1 <- as.matrix(simulate(m.full.ab, nsim=1000))
tmp2 <- aggregate(tmp1, by=list(spp=fern.data.ab$spp), sum, na.rm=TRUE)
tmp3 <- apply(tmp2[,-1], 2, sort, decreasing=TRUE)
tmp4 <- aggregate(fern.data.ab$abundance, by=list(spp=fern.data.ab$spp), sum)
plot(rad(tmp2[,2]), type="n")
for(i in sample(ncol(tmp2),100))
    lines(rad(tmp2[,i]), col="grey")
lines(rad(tmp4$x))

f1 <- function(.){
    pred <- predict(., type="response")
    aggregate(pred, by=list(spp=fern.data.ab$spp),sum)$x
}

hist(predict(m.full.ab, random.only=TRUE))
