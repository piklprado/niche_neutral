# Libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)

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

write.csv(prop_df, "results/supplementary_01.csv")

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

variable_new <- c("Conditional", "Fixed", "Random", "(1+(G|L))", "(1|SP:L)", "(1|SP:R)",
                  "(1|SP)", "(1+(G|SP))", "(1|L)")
names(variable_new) <- unique(effects_df$variable)
effects_df$variable_name <- str_replace_all(effects_df$variable, variable_new) %>%
  factor(.,
         levels = c("Conditional", "Random",
                    "Fixed", "(1+(G|L))", "(1+(G|SP))", "(1|SP)",
                    "(1|SP:R)", "(1|SP:L)", "(1|L)"))

effects_df$scenario <- factor(effects_df$scenario,
                               levels = c("Stochastic with right traits",
                                          "Deterministic with right traits",
                                          "Deterministic with wrong traits"))



# type of effect
type <- data.frame(variable = unique(effects_df$variable),
                   type = c("Total", "Selection", "Random", "Selection", "Drift", "Drift",
                            "Selection", "Selection", "Drift"))

df <- left_join(effects_df, type, by = "variable") %>%
  filter(!type %in% c("Random", "Total"))

p <- ggplot(data = effects_df, aes(y = variable, x = value, group = scenario)) +
  geom_point(aes(color = scenario)) +
  geom_segment(aes(y = variable, yend = variable,
                   x = lwr, xend = upr, color = scenario)) +
  # #facet_wrap(~ scenario, ncol = 1) +
  labs(x = "Partitioned R-squared value", y = "Model term") +
  theme_classic()

p <- ggplot(data = df, aes(y = variable_name, x = value, shape = scenario)) +
  #geom_jitter(aes(color = type)) +
  geom_point(aes(color = type)) +
  scale_shape_manual(values = c(19, 17, 2)) +
  geom_segment(aes(y = variable_name, yend = variable_name,
                   x = lwr, xend = upr, color = type)) +
  scale_color_manual(values = c(neutral, niche)) +
  labs(x = "Partitioned R-squared value", y = "Model term") +
  theme_classic()


p

png("figures/S1.png", res = 300, width = 1800, height = 1200)
p
dev.off()
