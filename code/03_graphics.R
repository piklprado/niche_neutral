#####################################################################################
# Partitioning niche and neutral dynamics on community assembly #####################
# Mortara et al #####################################################################
#####################################################################################

## Code for making graphics
## Sara Mortara and Paulo Inacio Prado

## loading packages
library(dplyr)
library(ggplot2)
library(gridExtra)

## functions to make graphics
source("functions.R")


## reading data
all.preds <- read.csv("../results/predicted.csv")
sads.meta <- read.csv("../results/sads_metacommunity.csv")
fern.data <- read.csv("../data/data_for_modeling.csv")

head(sads.meta)

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
    fig.meta1 + scale_colour_Publication() + theme_Publication() +,
    fig.meta2 + scale_colour_Publication() + theme_Publication() + theme(legend.position = "none"),
    fig.meta3 + scale_colour_Publication() + theme_Publication()
)

## Arrange with grid.arrange
cairo_pdf("../figures/rad_mettacomunity.pdf", width = 9, height = 7.5)

grid.arrange(
    grobs=fig.meta.list,
    bottom=textGrob("Altitude (m)", gp=gpar(fontface = "bold", cex=1.2), vjust=-1.25),
    layout_matrix=matrix(c(1,1,2,3), ncol=2, byrow=TRUE)
)

dev.off()
