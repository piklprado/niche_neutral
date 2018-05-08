######################################################################################
# Tailoring species abundance distributions with trait-environmental correlations ####
# Mortara et al ######################################################################
#####################################################################################

# Code for applying the niche-neutral GLMM framework

#############################
# PART 1: loading packages #
############################

# required packages
library(bbmle) 
library(lme4)
library(optimx)


#############################
# PART 2: loading data ######
############################

fern.data <- read.csv("fern_data_Mortaraetal.csv", header=TRUE)
head(fern.data)

#####################################################################
# PART 3: building models to represent hypothesis ##################
# Step by step building models corresponding to general hypothesis #
####################################################################

##################################
# 1. Idiosyncratic hypothesis ###
#################################
# 1.1 Null model
m1 <- glmer(abundance ~ 1 + (1|mountain) + (1|species) , data=fern.data,
                           family="poisson")

###################################
# 2. Neutral dynamic hypothesis ##
###################################

# 2.1 Neutral model with species restricted to mountains 
m2.1 <- glmer(abundance ~ 1 + (1|mountain:species), data=fern.data, family="poisson")

# 2.2 Neutral model with species restricted to mountains and sampling sites
m2.2 <- glmer(abundance ~ 1 + (1|mountain:species)+ (1|alt_std:species), data=fern.data, family="poisson") 

########################
# 3. Niche hypothesis #
#######################

##########################################
# 3.1. Absolute fitness among strategies #
##########################################

# 3.1.1 Ecological Strategy defined as combination of laminar tickness and indumentum
m3.1.1 <- glm(abundance ~ thickness + indumentum, data=fern.data, family="poisson")

# 3.1.2 Ecological Strategy defined as a combination of indumendum and life form
m3.1.2 <- glm(abundance ~ indumentum  + life_form,
            data=fern.data, family="poisson")

# 3.1.3. Ecological Strategy defined as a combination of laminar tickness and life form
m3.1.3 <- glm(abundance ~ thickness 
            +  life_form,
            data=fern.data, family="poisson")

# 3.1.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum
m3.1.4 <- glm(abundance ~ thickness + indumentum
            +  life_form,
            data=fern.data, family="poisson")

########################################
# 3.2. Trait-environment correlation ##
########################################

# 3.2.1 Ecological Strategy defined as combination of laminar tickness and indumentum interacting with altitude
m3.2.1 <- glm(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2),
              data=fern.data, family="poisson")

# 3.2.2 Ecological Strategy defined as a combination of indumendum and life form interacting with altitude
m3.2.2 <- glm(abundance ~ indumentum*alt_std + indumentum*I(alt_std^2)
              + life_form*alt_std + life_form*I(alt_std^2),
            data=fern.data, family="poisson")

# 3.2.3. Ecological Strategy defined as a combination of laminar tickness and life form interacting with altitude
m3.2.3 <- glm(abundance ~ thickness*alt_std + thickness*I(alt_std^2) 
            +  life_form*alt_std + life_form*I(alt_std^2),
            data=fern.data, family="poisson")

# 3.2.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude
m3.2.4 <- glm(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2),
            data=fern.data, family="poisson")


################################
# 4. Emergent Group hypothesis #
################################

###################################
# 4.1 Absolute fitness with drift #
###################################

# 4.1.1 Ecological Strategy defined as combination of laminar tickness and indumentum, drift among species sharing the same ES
m4.1.1 <- glmer(abundance ~ thickness + indumentum + (1|species), data=fern.data, family="poisson")

# 4.1.2 Ecological Strategy defined as a combination of indumendum and life form, drift among species sharing the same ES
m4.1.2 <- glmer(abundance ~ indumentum  + life_form + (1|species),
            data=fern.data, family="poisson")

# 4.1.3. Ecological Strategy defined as a combination of laminar tickness and life form, drift among species sharing the same ES
m4.1.3 <- glmer(abundance ~ thickness 
            +  life_form + (1|species),
            data=fern.data, family="poisson")

# 4.1.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum, drift among species sharing the same ES
m4.1.4 <- glmer(abundance ~ thickness + indumentum
            +  life_form + (1|species),
            data=fern.data, family="poisson")

##################################################
# 4.2. Trait-environment correlation with drift ##
##################################################

# 4.2.1 Ecological Strategy defined as combination of laminar tickness and indumentum interacting with altitude, drift among species sharing the same ES
m4.2.1 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
             + (1|species),
              data=fern.data, family="poisson",
             nAGQ=10, control=glmerControl(optimizer="bobyqa"))

# 4.2.2 Ecological Strategy defined as a combination of indumendum and life form interacting with altitude, drift among species sharing the same ES
m4.2.2 <- glmer(abundance ~ indumentum*alt_std + indumentum*I(alt_std^2)
              + life_form*alt_std + life_form*I(alt_std^2)
             + (1|species),
            data=fern.data, family="poisson",
            nAGQ=10, control=glmerControl(optimizer="bobyqa"))

# 4.2.3. Ecological Strategy defined as a combination of laminar tickness and life form interacting with altitude, drift among species sharing the same ES
m4.2.3 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2) 
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species),
            data=fern.data, family="poisson",
           nAGQ=10, control=glmerControl(optimizer="bobyqa"))

# 4.2.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES
m4.2.4 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species),
            data=fern.data, family="poisson",
            nAGQ=10, control=glmerControl(optimizer="bobyqa"))

#######################################################
# 5. Emergent Group with limited dispersal hypothesis #
#######################################################

##################################################################
# 5.1 Absolute fitness with drift and regional limited dispersal #
##################################################################

# 5.1.1 Ecological Strategy defined as combination of laminar tickness and indumentum, drift among species sharing the same ES, regional limited dispersal
m5.1.1 <- glmer(abundance ~ thickness + indumentum + (1|species) + (1|mountain:species), data=fern.data, family="poisson")

# 5.1.2 Ecological Strategy defined as a combination of indumendum and life form, drift among species sharing the same ES, regional limited dispersal
m5.1.2 <- glmer(abundance ~ indumentum  + life_form + (1|species)+ (1|mountain:species),
            data=fern.data, family="poisson")

# 5.1.3. Ecological Strategy defined as a combination of laminar tickness and life form, drift among species sharing the same ES, regional limited dispersal
m5.1.3 <- glmer(abundance ~ thickness 
            +  life_form + (1|species) + (1|mountain:species),
            data=fern.data, family="poisson")

# 5.1.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum, drift among species sharing the same ES, regional limited dispersal
m5.1.4 <- glmer(abundance ~ thickness + indumentum
            +  life_form + (1|species) + (1|mountain:species),
            data=fern.data, family="poisson")

##########################################################################
# 5.2 Absolute fitness with drift, local  and regional limited dispersal #
##########################################################################

# 5.2.1 Ecological Strategy defined as combination of laminar tickness and indumentum, drift among species sharing the same ES, local and regional limited dispersal
m5.2.1 <- glmer(abundance ~ thickness + indumentum + (1|species) +
                    (1|mountain:species) + (1|alt_std:species),
              data=fern.data, family="poisson")

# 5.2.2 Ecological Strategy defined as a combination of indumendum and life form, drift among species sharing the same ES, local and regional limited dispersal
m5.2.2 <- glmer(abundance ~ indumentum  + life_form + (1|species)+
                    (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson")

# 5.2.3. Ecological Strategy defined as a combination of laminar tickness and life form, drift among species sharing the same ES, local and regional limited dispersal
m5.2.3 <- glmer(abundance ~ thickness 
            +  life_form + (1|species) +
            (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson")

# 5.2.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum, drift among species sharing the same ES, local and regional limited dispersal
m5.2.4 <- glmer(abundance ~ thickness + indumentum
            +  life_form + (1|species) +
            (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson")


################################################################################
# 5.3 Trait-environment correlation with drift and regional limited dispersal ##
################################################################################

# 5.3.1 Ecological Strategy defined as combination of laminar tickness and indumentum interacting with altitude, drift among species sharing the same ES, regional limited dispersal
m5.3.1 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
             + (1|species) + (1|mountain:species),
              data=fern.data, family="poisson",
             nAGQ=1, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))

# 5.3.2 Ecological Strategy defined as a combination of indumendum and life form interacting with altitude, drift among species sharing the same ES, regional limited dispersal
m5.3.2 <- glmer(abundance ~ indumentum*alt_std + indumentum*I(alt_std^2)
              + life_form*alt_std + life_form*I(alt_std^2)
             + (1|species) + (1|mountain:species),
            data=fern.data, family="poisson",
            nAGQ=1, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
            

# 5.3.3. Ecological Strategy defined as a combination of laminar tickness and life form interacting with altitude, drift among species sharing the same ES, regional limited dispersal
m5.3.3 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2) 
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1|mountain:species),
            data=fern.data, family="poisson",
           nAGQ=1, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))

# 5.3.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, regional limited dispersal
m5.3.4 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1|mountain:species),
            data=fern.data, family="poisson",
           control=glmerControl(optimizer="optimx",
                optCtrl=list(method="nlminb")))

                                        #nAGQ=1, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)),
           #control=glmerControl(optimizer="optimx",
            #    optCtrl=list(method="nlminb")))
#nAGQ=1)

########################################################################################
# 5.4 Trait-environment correlation with drift, local and regional limited dispersal ##
#######################################################################################

# 5.4.1 Ecological Strategy defined as combination of laminar tickness and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m5.4.1 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
             + (1|species) + (1|mountain:species) + (1|alt_std:species),
              data=fern.data, family="poisson",
             control=glmerControl(optimizer="optimx",
                optCtrl=list(method="nlminb")))

# 5.4.2 Ecological Strategy defined as a combination of indumendum and life form interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m5.4.2 <- glmer(abundance ~ indumentum*alt_std + indumentum*I(alt_std^2)
              + life_form*alt_std + life_form*I(alt_std^2)
             + (1|species) + (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson",
            control=glmerControl(optimizer="optimx",
                optCtrl=list(method="nlminb")))

# 5.4.3. Ecological Strategy defined as a combination of laminar tickness and life form interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m5.4.3 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2) 
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson",
           control=glmerControl(optimizer="optimx",
                optCtrl=list(method="nlminb")))

# 5.4.4 Ecological Strategy defined by all the three traits: laminar thickness, life form and indumentum interacting with altitude, drift among species sharing the same ES, local and regional limited dispersal
m5.4.4 <- glmer(abundance ~ thickness*alt_std + thickness*I(alt_std^2)
              + indumentum*alt_std + indumentum*I(alt_std^2)
            +  life_form*alt_std + life_form*I(alt_std^2)
           + (1|species) + (1|mountain:species) + (1|alt_std:species),
            data=fern.data, family="poisson",
           control=glmerControl(optimizer="optimx",
                optCtrl=list(method="nlminb")))

#            nAGQ=1)#, control=glmerControl(optimizer="bobyqa"))


#####################################################################
# PART 4: Model selection                          ##################
####################################################################

aic.tab <- AICtab(m1, # Null
       m2.1, m2.2, # Neutral
       m3.1.1, m3.1.2, m3.1.3, m3.1.4,
       m3.2.1, m3.2.2, m3.2.3, m3.2.4, # Niche
       m4.1.1, m4.1.2, m4.1.3, m4.1.4, 
       m4.2.1, m4.2.2, m4.2.3, m4.2.4, # Emergent Groups
       m5.1.1, m5.1.2, m5.1.3, m5.1.4, 
       m5.2.1, m5.2.2, m5.2.3, m5.2.4, # EG, drift, regional limited dispersal
       m5.3.1, m5.3.2, m5.3.3, m5.3.4,
       m5.4.1, m5.4.2, m5.4.3, m5.4.4, # EG, drift, local and regional limited dispersal
       weights=TRUE, base=TRUE) 

aic.tab

library(xtable)
aic.tab <- print(aic.tab, min.weight=0.001)

print(xtable(aic.tab, caption="Model selection confronting idiosyncratic, ecological drift, niche, emergent groups,
and emergent groups with ecological drift models for fern species abundances on local communities."), caption.placement="top", label="aic")


#####################################################################
# PART 5: Calculating predicted values from best model  #############
####################################################################

# First, we create a data frame with all combination of sites, species and traits 
comb.table <- data.frame(expand.grid(mountain=levels(fern.data$mountain),
                        alt_std=unique(fern.data$alt_std),
                        species=unique(fern.data$species)),
                       life_form=fern.data$life_form,
                       thickness=fern.data$thickness)

comb.table <- na.omit(comb.table)

# Second, we use the function predict to create a data frame of predicted values for all possible combinations based on the best model m5.4.3
pred.values <- predict(m5.4.3, re.form=NULL, newdata=comb.table,
                       type='response')

# Third we calculate mean predicted values and standard error for each altitude 

## Predicted mean values
pred.table <- aggregate(pred.values, by=list(altitude=comb.table$alt_std, thickness=comb.table$thickness,
                                             life_form=comb.table$life_form), mean)

names(pred.table)[4] <-  "mean"

## Predicted stardard error
pred.table$se <- aggregate(pred.values, by=list(altitude=comb.table$alt_std, thickness=comb.table$thickness,
                                             life_form=comb.table$life_form),
                            function(x)sd(x)/sqrt(length(x)))$x

head(pred.table)

# Finally, we calculate the upper and lower confidence interval based on t distribution 
## Confidence Interval (mean +- standard error * t(pdf)
t.prev <- pt(pred.table$mean, df=(nrow(pred.table)-1))
pred.table$lower <- (pred.table$mean - pred.table$se)*t.prev
pred.table$upper <- (pred.table$mean + pred.table$se)*t.prev


#####################################################################
# PART 6: Calculating R-squared for the best model  ################
####################################################################
library(piecewiseSEM)

r2 <- sem.model.fits(m5.4.3)
r2

#marginal 0.03
#conditional 0.86

# Second we create a data frame with observed mean values and its standard error
obs <- aggregate(fern.data$abundance, by=list(altitude=fern.data$alt_std, thickness=fern.data$thickness,
                                              life_form=fern.data$life_form), mean)
## Observed standard error
obs$se <- aggregate(fern.data$abundance, by=list(altitude=fern.data$alt_std, thickness=fern.data$thickness,
                                              life_form=fern.data$life_form),
                            function(x)sd(x)/sqrt(length(x)))$x
names(obs) <- c("Altitude", "thickness", "life_form", "Mean Abundance", "std")

#save.image("Mortaraetal.RData")


    

