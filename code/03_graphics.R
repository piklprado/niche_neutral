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
library(ggpubr)

## functions to make graphics
source("code/functions.R")

## reading data
all.preds <- read.csv("results/predicted.csv")
sads.meta <- read.csv("results/sads_metacommunity.csv")
fern.data <- read.csv("data/data_for_modeling.csv")
spp <- read.csv("data/species_list.csv")

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
plot(1:3, col = paleta2, pch = 19, cex = 5)

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
fig.meta2 <- fern.data %>%
    filter(spp == "sp120") %>% # rank 1
    ggplot(aes(altitude, abundance)) +
    stat_smooth(se = FALSE, aes(color = region), alpha = 0.5, span = 0.6) +
    geom_jitter(aes(color = region), size = 2.5, alpha = 0.5) +
    scale_color_manual(values = paleta2) +
    labs(title = expression(bold("Core:") ~ italic("Polybotrya cylindrica")),
                            x = "", y = "Abundance", color = "Region", tag = "B") +
    theme_classic()

#fig.meta2

## 3rd panel: Most abundant of occasional species
fig.meta3 <- fig.meta2 %+% filter(fern.data, spp == "sp159") +
    labs(title = expression(bold("Occasional:") ~ italic("Trichomanes polypodioides")),
         x = "", y = "", color = "Region", tag = "D") +
    theme_classic()

#fig.meta3

# New proposal, selection of six species to show differences among life forms
# selected species (ep, hemi, ter)
# Occasional: sp159, sp033, sp028
# Rare: sp084, sp120, sp012
# Core ep

sp.sel <- c("sp158", "sp033", "sp028",
            "sp084", "sp120", "sp030")


names.sp <- spp %>%
    filter(species %in% sp.sel) %>%
    select(name_accepted, species) %>%
    rename(spp = species)

fig.sp <- list()

data.sel <- fern.data %>%
    filter(spp %in% sp.sel) %>%
    left_join(names.sp, by = "spp") %>%
    mutate(spp = factor(.$spp,
                        levels = c("sp084", "sp120", "sp030",
                                   "sp158", "sp033", "sp028"))) %>%
    arrange(spp)


# fig.sp <- data.sel %>%
#     ggplot(aes(altitude, abundance)) +
#     stat_smooth(se = FALSE, aes(color = region), alpha = 0.5, span = 0.6) +
#     geom_jitter(aes(color = region), size = 2.5, alpha = 0.5) +
#     scale_color_manual(values = paleta2) +
#     facet_wrap(~ name_accepted, scales = "free",
#                nrow = 3, ncol = 2) +
#     labs(#title = expression(bold("Core:") ~ italic("Polybotrya cylindrica")),
#         x = "Altitude (m)", y = "Abundance", color = "Region") +
#     theme_classic() +
#     theme(strip.text = element_text(hjust = 0, face = "italic"),
#           strip.background = element_blank(),
#           legend.position = "bottom")

fig.sp <- list()

for (i in 1:length(sp.sel)) {
    fig.sp[[i]] <- ggplot(data = filter(data.sel, spp == unique(data.sel$spp)[i]),
                          aes(altitude, abundance)) +
        stat_smooth(se = FALSE, aes(color = region), alpha = 0.5, span = 0.6) +
        geom_jitter(aes(color = region), size = 2.5, alpha = 0.5) +
        scale_color_manual(values = paleta2) +
        labs(title = "",
             subtitle = unique(data.sel$name_accepted)[i],
             x = "", y = "", color = "Region",
             tag = LETTERS[i]) +
        theme_classic() +
        theme(plot.subtitle = element_text(face = "italic"))
}


# writing text around plots
fig1 <- annotate_figure(fig.sp[[1]] + theme(legend.position = "none"),
                        right = text_grob(bquote(""), rot = -90))
fig2 <- annotate_figure(fig.sp[[2]] + theme(legend.position = "none"),
                        right = text_grob(bquote(""), rot = -90))
fig3 <- annotate_figure(fig.sp[[3]] + theme(legend.position = "none"),
                        right = text_grob(bquote(bold("Core")), rot = -90))
fig6 <- annotate_figure(fig.sp[[6]] + theme(legend.position = "none"),
                        right = text_grob(bquote(bold("Occasional")), rot = -90))
# fig6 <- annotate_figure(fig.sp[[6]] + theme(legend.position = "none"),
#                         right = text_grob(bquote(bold("Terrestrial")), rot = -90))

# combining species figures
fig.sp.combined <- ggarrange(fig1 +
                                 ggtitle(expression(bold("Epiphyte"))) +
                                 theme(plot.title = element_text(hjust = 0.5, size = 12)),
                             fig2 +
                                 ggtitle(expression(bold("Hemiepiphyte"))) +
                                 theme(plot.title = element_text(hjust = 0.5, size = 12)),
                             fig3 +
                                 ggtitle(expression(bold("Terrestrial"))) +
                                 theme(plot.title = element_text(hjust = 0.5, size = 12)),
                             fig.sp[[4]],
                             fig.sp[[5]],
                             fig6,
                             ncol = 3, nrow = 2,
                             common.legend = TRUE, legend = "bottom")

annotate_figure(fig.sp.combined, left = "Abundance", bottom = "Altitude")

# saving combined plot of species
cairo_pdf("figures/core_ocasional_species.pdf", width = 8.5, height = 7.5)

annotate_figure(fig.sp.combined, left = "Abundance", bottom = "Altitude")

dev.off()


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

fig.meta1 + #scale_colour_Publication() +
    labs(tag = "") +
    theme_Publication()
# grid.arrange(
#     grobs = fig.meta.list,
#     bottom = textGrob("Altitude (m)", gp = gpar(fontface = "bold", cex = 1.2), vjust = -1.25),
#     layout_matrix = matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE)
# )

dev.off()

# Figure R2 ####
# it wont be a figure, just values on text

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
