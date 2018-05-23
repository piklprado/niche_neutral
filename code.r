######################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al ######################################################################
#####################################################################################

## Code for applying the niche-neutral GLMM framework
## Data and analysis from the manuscript  

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
source("r2_table.R")

#############################
# PART 2: loading data ######
############################

fern.data <- read.csv("fern_data_Mortaraetal.csv", header=TRUE)
head(fern.data)

fern.data$site <- scale(rep(1:30, length(unique(fern.data$species))))


#####################################################################
# PART 3: building the model to represent our hypothesis ############
# Step by step building models corresponding to general hypothesis #
####################################################################

# Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m.full <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              #+ indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1|species:mountain) + (1|species:site) + (1 + alt_std+I(alt_std^2)|species),
           data=fern.data, family="poisson",
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.neutral <- glmer(abundance ~ (1|species) + (1|species:mountain) + (1|species:site) + (1 + alt_std|species),
           data=fern.data, family="poisson",
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e6)))

m.niche <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              #+ indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1 + alt_std|species),
           data=fern.data, family="poisson",
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.null <- glmer(abundance ~ 1 +
           (1|species) + (1|mountain) + (1|site) + (1 + alt_std|species),
           data=fern.data, family="poisson",
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e6)))

m.list <- list(m.full, m.neutral, m.niche, m.null)

bic.tab <- sapply(m.list, BIC)
mod.names <- c("niche & neutral", "neutral", "niche", "null") 
names(bic.tab) <- mod.names

bic.tab

r2.tab <- sapply(m.list, r2.table)

r2.tab <- bind_rows(r2.tab)

row.names(r2.tab) <- mod.names

r2.tab

#####################################################################
# PART 4: Calculating predicted values from best model  #############
####################################################################

# First, we create a data frame with all combination of sites, species and traits 
comb.table <- data.frame(expand.grid(mountain=levels(fern.data$mountain),
                                     alt_std=unique(fern.data$alt_std),
                                     site=unique(fern.data$site),
                                     species=unique(fern.data$species)),
                         life_form=fern.data$life_form,
                         thickness=fern.data$thickness,
                         indumentum=fern.data$indumentum)

comb.table <- na.omit(comb.table)

# Second, we use the function predict to create a data frame of predicted values for all possible combinations based on the best model m5.4.3
pred.values <- predict(m.full, re.form=NULL, newdata=comb.table,
                       type='response')

# Third we calculate mean predicted values and standard error for each altitude 

## Predicted mean values
pred.table <- aggregate(pred.values, list(altitude=comb.table$alt_std, thickness=comb.table$thickness,
                                          indumentum=comb.table$indumentum,
                                             life_form=comb.table$life_form), mean)

names(pred.table)[5] <- "mean"


## Predicted stardard error
pred.table$se <- aggregate(pred.values, by=list(altitude=comb.table$alt_std, thickness=comb.table$thickness,
                                                indumentum=comb.table$indumentum,
                                             life_form=comb.table$life_form),
                            function(x)sd(x)/sqrt(length(x)))$x

head(pred.table)

# Finally, we calculate the upper and lower confidence interval based on t distribution 
## Confidence Interval (mean +- standard error * t(pdf)
t.prev <- pt(pred.table$mean, df=(nrow(pred.table)-1))
pred.table$lower <- (pred.table$mean - pred.table$se)*t.prev
pred.table$upper <- (pred.table$mean + pred.table$se)*t.prev

# Second we create a data frame with observed mean values and its standard error
obs <- aggregate(fern.data$abundance, by=list(altitude=fern.data$alt_std, thickness=fern.data$thickness,
                                              indumentum=fern.data$indumentum,
                                              life_form=fern.data$life_form), mean)
## Observed standard error
obs$se <- aggregate(fern.data$abundance, by=list(altitude=fern.data$alt_std, thickness=fern.data$thickness,
                                                 indumentum=fern.data$indumentum,
                                              life_form=fern.data$life_form),
                            function(x)sd(x)/sqrt(length(x)))$x

head(obs)
names(obs) <- c("Altitude", "thickness", "indumentum", "life_form", "Abundance", "std")

#############################################
######### Creating figures ##################
############################################

###############################################
######### DAQUI PRA FRENTE AINDA NAO FUNFA ####
###############################################

# all trait combinations
esp.hab <- expand.grid(c('membranacea', 'coriacea'), c('ausente','presente'), c('ter', 'hemi', 'ep'))
esp.hab

#########################################
#### GRAFICO MODELO #####################
#########################################

cor1 <-rgb(140, 1, 28, maxColorValue=255) #rgb(44, 152, 32, maxColorValue=255) # terrestre
cor3 <- rgb(4, 70, 120, maxColorValue=255) #rgb(239, 144, 33, maxColorValue=255) # hemi
cor2 <- rgb(199, 172, 29, maxColorValue=255) # ep

head(obs)

ep.cor.si <- subset(obs, thickness=="coriacea" & indumentum=="ausente" &life_form=="ep")
ep.cor.ci <- subset(obs, thickness=="coriacea" & indumentum=="presente" &life_form=="ep")
ep.mem.si <- subset(obs, thickness=="membranacea" & indumentum=="ausente" &life_form=="ep")
ep.mem.ci <- subset(obs, thickness=="membranacea" & indumentum=="presente" &life_form=="ep")

head(obs)

par(mfrow=c(1,2))
plot(Abundance ~ Altitude, data=ep.cor.si, log='x', ylim=c(0,20), col=cor1, las=1)
segments(x0=ep.cor.si[,1],
         y0= ep.cor.si[,5] + ep.cor.si[,6],
         y1= ep.cor.si[,5] - ep.cor.si[,6], col=cor1)
points(Abundance ~ Altitude, data=ep.cor.ci, pch=19,col=cor1)
segments(x0=ep.cor.ci[,1],
         y0= ep.cor.ci[,5] + ep.cor.ci[,6],
         y1= ep.cor.ci[,5] - ep.cor.ci[,6], col=cor1)
plot(Abundance ~ Altitude, data=ep.mem.si, log='x', ylim=c(0,20), col=cor1, las=1)
segments(x0=ep.mem.si[,1],
         y0= ep.mem.si[,5] + ep.mem.si[,6],
         y1= ep.mem.si[,5] - ep.mem.si[,6], col=cor1)
points(Abundance ~ Altitude, data=ep.mem.ci, pch=19, col=cor1)
segments(x0=ep.mem.ci[,1],
         y0= ep.mem.ci[,5] + ep.mem.ci[,6],
         y1= ep.mem.ci[,5] - ep.mem.ci[,6], col=cor1)


head(ep.cor.ci)


loadfonts(device = "postscript")

pdf("graf_modelo.pdf")

par(mai=c(0.5, 0.5, 0.2, 0.25), oma=c(1, 1, 1, 0.1))
layout(matrix(c(0, 0, 0, 0,
                0, 1, 2, 0,
                0, 3, 4, 0,
                0, 5, 6, 0),4,4, byrow=TRUE),
                widths=c(0.1, 1, 1, 0.1), heights=0.1)

for(i in 1:8){
plot(obs[obs$thickness==esp.hab[i,1] & obs$indumentum==esp.hab[i,2] & obs$life_form==esp.hab[i,3], c(1,5)],
     pch=rep(c(21, 19), 3)[i], bty="l", xlab="", ylab="", cex=1.7, yaxt="n", xaxt="n", log='y')
     #pt.bg='white' )#, 
     col=rep(c(cor1, cor2, cor3), each=2)[i], 
     ylim=rbind(c(0.1,40), c(0.1,40), c(0.1,120), c(0.1,120), c(0.1,23), c(0.1,23))[i,])
# controlando eixos
if(i %in% c(5,6)){
    axis(1, at=unique(com.obs$Altitude), labels=unique(com.obs$Altitude), cex.axis=1.3)}
else{axis(1, at=unique(com.obs$Altitude), labels=FALSE)}
if(i %in% seq(1,6,2)){
    axis(2, cex.axis=1.3, las=1)}
else{axis(2, labels=FALSE)}
# erro padrao obs
segments(x0=com.obs[com.obs$esp==esp.hab[i,1]  & com.obs$hab==esp.hab[i,2], 1],
         y0=com.obs[com.obs$esp==esp.hab[i,1]  & com.obs$hab==esp.hab[i,2], 4] +
         com.obs[com.obs$esp==esp.hab[i,1]  & com.obs$hab==esp.hab[i,2], 5],
         y1=com.obs[com.obs$esp==esp.hab[i,1]  & com.obs$hab==esp.hab[i,2], 4] -
         com.obs[com.obs$esp==esp.hab[i,1]  & com.obs$hab==esp.hab[i,2], 5],
         col=rep(c(cor1, cor2, cor3), each=2)[i])
## Previsto medio
lines(com.prev[com.prev$esp==esp.hab[i,1] & com.prev$hab==esp.hab[i,2], c(1,4)],
      col=rep(c(cor1, cor2, cor3), each=2)[i])
## Intervalo de mais ou mesno 2 x se
lines(com.prev[com.prev$esp==esp.hab[i,1] &  com.prev$hab==esp.hab[i,2], c(1,6)], lty=2,
     col=rep(c(cor1, cor2, cor3), each=2)[i])
lines(com.prev[com.prev$esp==esp.hab[i,1] &  com.prev$hab==esp.hab[i,2], c(1,7)], lty=2,
      col=rep(c(cor1, cor2, cor3), each=2)[i])
mtext(paste(paste("(", letters[1:6][i], sep=""), ")", sep=""), side=3, adj=0.05, padj=-0.5, cex=1) #font=2
}
mtext("Mean species abundances (log)", side=2, outer=TRUE, padj=1, cex=1.2)
mtext("Altitude (m)", side=1, outer=TRUE, padj=-0.5, cex=1.2)
mtext("Membranaceous", side=3, adj=0.25, padj=1, outer=TRUE, font=2)
mtext("Coriaceous", side=3, outer=TRUE, adj=0.8, padj=1, font=2)

mtext("Terrestrial", side=4, outer=TRUE, padj=-1.7, adj=0.87, font=2)
mtext("Hemiepiphyte", side=4, outer=TRUE, padj=-1.7, font=2)
mtext("Epiphyte", side=4, outer=TRUE, padj=-1.7, adj=0.135, font=2)

dev.off()

embed_fonts("graf_modelo.eps", outfile = "graf_modelo.eps",
            options = "-dEPSCrop")


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


    

