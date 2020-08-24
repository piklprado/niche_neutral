# Libraries
library(ggplot2)
library(dplyr)
library(reshape2)

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
models <- c("Drift-Selection", "Drift-Selection (no traits)", "Selection", "Selection (no traits)", "Drift", "Idyosyncratic")

# Proportions of the best fitted models
prop_list <- lapply(prop_files, read.csv)
prop_df <- t(bind_cols(lapply(prop_list, function(x) x[, 2])))
colnames(prop_df) <- models
rownames(prop_df) <- comm

write.table(prop_df, "results/supplementary_01.csv")

# Effects -----------------------------------------------------------------
effects_list <- lapply(effects_files, read.csv)
names(effects_list) <- comm
effects_all <- bind_rows(effects_list, .id = "scenario") %>%
  rename(stat = X)
effects_all$stat[effects_all$stat == ""] <- "mean"

effects_mean <- melt(effects_all) %>%  filter(stat == "mean")
lwr <- melt(effects_all) %>%  filter(stat == "2.5%") %>% select(value) %>% rename(lwr = value)
upr <- melt(effects_all) %>%  filter(stat == "97.5%") %>% select(value) %>%
  rename(upr = value)
effects_df <- cbind(effects_mean, lwr, upr)
effects_df$variable <- gsub("X.1.", "", effects_df$variable, fixed = TRUE)
effects_df$scenario <- factor(effects_df$scenario,
                               levels = c("Stochastic with right traits",
                                          "Deterministic with right traits",
                                          "Deterministic with wrong traits"))

p <- ggplot(data = effects_df, aes(y = variable, x = value, group = scenario)) +
  geom_point(aes(color = scenario)) +
  geom_segment(aes(y = variable, yend = variable,
                   x = lwr, xend = upr, color = scenario)) +
  # #facet_wrap(~ scenario, ncol = 1) +
  labs(x = "Partitioned R-squared value", y = "Model term") +
  theme_classic()

p

png("figures/S1.png", res = 300, width = 1800, height = 1200)
p
dev.off()

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
