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
source("code/functions.R")

# Colors from manuscript
# niche <- rgb(202, 0, 32, maxColorValue = 255)
# neutral <- rgb(4, 4, 160, maxColorValue = 255)
# nineu <- rgb(245, 219, 50, maxColorValue = 255)
grey <- rgb(180, 177, 172, maxColorValue = 255)

neutral <- "#091346"
niche <- "#801411"
nineu <- "#D6A838" #"#e1b727"

# trying colors for mountains, but ended with kind of grey scale
GR <- '#000000'
MR <- '#2a623d' #or '#5d5d5d'
PR <- '#aaaaaa'

paleta <- c(niche, neutral, nineu, grey)
paleta2 <- c(GR, MR, PR)

# plot(1:4, col = paleta, pch = 19, cex = 5)
# plot(1:4, col = c(nineu, paleta2), pch = 19, cex = 5)

## reading data
all.preds <- read.csv("results/predicted.csv")
sads.meta <- read.csv("results/sads_metacommunity.csv")
fern.data <- read.csv("data/data_for_modeling.csv")

head(sads.meta)

### The figure ###
## 1st panel of the figure: RAD with total and random part of standard deviations
fig.meta1 <- sads.meta %>%
    ggplot(aes(x = sp.rank, y = mean.obs)) +
    geom_ribbon(aes(ymin = lwr.obs, ymax = upr.obs), alpha = 0.3, fill = grey)  +
    geom_linerange(aes(x = sp.rank, ymin = lwr.random, ymax = upr.random, color = ab.class), size = 0.4) +
    geom_point(aes(color = ab.class)) +
    scale_color_manual(values = c(nineu, neutral)) +
    labs(x = "Abundance rank", y = "Mean abundance", color = "", tag = "A") +
    scale_y_log10() +
    theme_classic()

fig.meta1

## 2nd panel: Abundances of selected species over the Gradient
## Most abundant of core species
fig.meta2 <- fern.data[fern.data$spp == sads.meta$spp[sads.meta$sp.rank == 1], ] %>%
    ggplot(aes(altitude, abundance)) +
    stat_smooth(se = FALSE, aes(color = region), alpha = 0.5, span = 0.6) +
    geom_jitter(aes(color = region), size = 2.5, alpha = 0.5) +
    scale_color_manual(values = paleta2) +
    labs(title = expression(bold("Core:") ~ italic("Polybotrya cylindrica")),
                            x = "", y = "Abundance", color = "Region", tag = "B") +
    theme_classic()

fig.meta2

## 3rd panel: Most abundant of occasional species
fig.meta3 <- fig.meta2 %+% fern.data[fern.data$spp == sads.meta$spp[sads.meta$sp.rank == 41],] +
    labs(title = expression(bold("Occasional:") ~ italic("Trichomanes polypodioides")),
                            x = "", y = "", color = "Region", tag = "C") +
    theme_classic()

fig.meta3

## Arranging all panels in a sigle figure, publication quality (see functions.R)
## List of panels
fig.meta.list <- list(
    fig.meta1 + #scale_colour_Publication() +
        theme_Publication(),
    fig.meta2 + #scale_colour_Publication() +
        theme_Publication() + theme(legend.position = "none"),
    fig.meta3 + #scale_colour_Publication() +
        theme_Publication()
)

## Arrange with grid.arrange
cairo_pdf("figures/rad_metacommunity.pdf", width = 9, height = 7.5)

grid.arrange(
    grobs = fig.meta.list,
    bottom = textGrob("Altitude (m)", gp = gpar(fontface = "bold", cex = 1.2), vjust = -1.25),
    layout_matrix = matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE)
)

dev.off()

# Figure R2

r2 <- read.csv("results/r2_table.csv")

r2$type[r2$type == "idiosyncratic"] <- "niche"

r2.df <- r2 %>%
    filter(!type %in% c("all", "all.random")) %>%
    mutate(R2 = round(R2, 4)) %>%
    group_by(group, type) %>%
    summarise(R2 = sum(R2))

#pdf("../figures/barplot.pdf")
ggplot(r2.df, aes(fill = type, y = R2, x = group)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic()
#dev.off()
