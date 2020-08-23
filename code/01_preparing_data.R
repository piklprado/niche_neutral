#####################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al #####################################################################
#####################################################################################

## Code for preparing data for modeling
### Separating abundant and rare spp

###############################
# PART 1: loading packages ####
###############################

library(reshape)

############################
# PART 2: loading data #####
############################

# reading data
fern.data <- read.csv("data/data_original.csv", as.is = TRUE)
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

sp.site <- apply(sp.site[, c(-1, -2)], 2, as.numeric)

dim(sp.site)

sp.alt <- apply(sp.site, 2, function(x) tapply(x, altitude, mean))

## Creating vector of species abundances for the total metacommunity
ab.meta <- colSums(sp.site)
names(ab.meta) <- sp

length(ab.meta)

ab.meta

#### separating 40 most abundant species in the metacommunity
indice.sp <- order(ab.meta, decreasing = TRUE)[1:40]

freq <- rowSums(sp.site > 0)

ab.mean <- apply(sp.site, 1, mean)
ab.mean

## Heat map of the species x sites matrix
image(x = 1:nrow(sp.site), y = 1:ncol(sp.site), log(sp.site[, order(-ab.meta)]))

#### Will indentify species as abundant (first 40 in rank) and rare (other)
ab.rare <- data.frame(spp = unique(fern.data$spp),
                      ab.rare = ifelse(sp %in% sp[indice.sp], "abundant", "rare"))

ab.rare

fern.data.new <- merge(fern.data, ab.rare)
dim(fern.data)
dim(fern.data.new)

head(fern.data.new)

fern.data.ab <- fern.data.new[fern.data.new$ab.rare == "abundant", ]
fern.data.rare <- fern.data.new[fern.data.new$ab.rare == "rare", ]

write.table(fern.data.new,
            "../data/data_for_modeling.csv",
            col.names = TRUE, row.names = FALSE, sep = ",")

write.table(fern.data.ab,
            "../data/data_ab_for_modeling.csv",
            col.names = TRUE, row.names = FALSE, sep = ",")

write.table(fern.data.rare,
            "../data/data_rare_for_modeling.csv",
            col.names = TRUE, row.names = FALSE, sep = ",")


