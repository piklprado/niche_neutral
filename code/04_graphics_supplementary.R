# Libraries
library(ggplot2)
library(dplyr)

#library(colortools)
# Setting color palette
neutral <- "#091346"
niche <- "#801411"
nineu <- "#D6A838" #"#e1b727"
idiosyncratic <- "#B4B1AC"

paleta <- c(niche, neutral, nineu, grey)

# Reading output from simulations
## Listing files
prop_files <- list.files("results/simulations", pattern = ".*_prop.csv$", full.names = TRUE)
effects_files <-  list.files("results/simulations", pattern = ".*_effects.csv$", full.names = TRUE)

# Proportion table --------------------------------------------------------
## Setting names
comm <- c("Deterministic with right traits", "Deterministic with wrong traits", "Stochastic with right traits")
models <- c("Niche-neutral", "Niche-neutral (no traits)", "Niche", "Niche (no traits)", "Neutral", "Idyosyncratic")

# Proportions of the best fitted models
prop_list <- lapply(prop_files, read.csv)
prop_df <- t(bind_cols(lapply(prop_list, function(x) x[, 2])))
colnames(prop_df) <- models
rownames(prop_df) <- comm


# Effects -----------------------------------------------------------------

effects_list <- lapply(effects_files, read.csv)
names(effects_list) <- comm
effects_df <- bind_rows(effects_list, .id = "scenario") %>%
  rename(stat = X)
effects_df$stat[effects_df$stat == ""] <- "mean"

effect_mean <- filter(effects_df, stat == "mean")
effect_lwr <- filter(effects_df, stat == "2.5%")
effect_lwr <- filter(effects_df, stat == "97.5%")



# other stuff w/ colors
# library(colortools)
# niche.cols <- analogous(niche)
# neutral.cols <- analogous(neutral)
# wheel(niche, num = 10)
# wheel(neutral, num = 10)
# wheel(nineu, num = 10)

# p1 <- ggplot(t1, aes(x = Community,
#                      y = Proportion,
#                      fill = Effect)) +
#   geom_bar(stat = "identity", color = "black") +
#   theme_classic()
#
#
# p1
# p2 <- ggplot(t2, aes(x = Community,
#                      y = Proportion,
#                      fill = Effect)) +
#   geom_bar(stat = "identity", color = "black") +
#   theme_classic()
#
# p2
#
# grid.arrange(p1, p2, nrow = 1)
#
# png("figures/test1.png")
# p1
# dev.off()
#
# png("figures/test2.png")
# p2
# dev.off()
#
# p1
# p2
